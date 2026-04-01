# SMZWWA_fixed.py
# ArcGIS Pro / ArcPy 3.x (Python 3.9+)
# Western WA workflow:
# 1) Type S & F buffers: (half bankfull + 50 ft), centerline-based
# 2) Type Np: ONLY Np connected to S/F; take FIRST HALF away from confluence toward upstream; buffer (half bankfull + 50 ft)
# 3) Type Np intersections (Western WA): fixed 56-ft buffers (true intersections only: >=2 distinct Np streams and >=1 interior contact)
# 4) Optional export of final merged buffers to shapefile (safe field mapping, optional RemoveZ/RemoveM, spatial index)
import arcpy
import os
import time

# --- Environment ---
arcpy.env.overwriteOutput = True
arcpy.env.addOutputsToMap = False  # avoid TOC/locks in Pro

# ---------------------------------------------------------------------
# User-provided parameters
# ---------------------------------------------------------------------
flowlines_shp = r"D:\Downloaded Shapes\NHD\Washington\DNR_Hydrography_-_Watercourses_SFNP_WWA_SJ.shp"
type_field = "FP_REF_ID"
bankfull_field = "BANKFULL_X"  # full bankfull width in feet
target_sr_wkid = 2913  # WA North (ftUS); use 2914 for WA South
snap_tol_ft = 10.0
intersection_radius_ft = 56.0  # fixed Western WA site radius
is_half_width = False
out_gdb = r"D:\Downloaded Shapes\NHD\Washington\WABuffers.gdb"

# LIKE patterns
s_like = "%S%"
f_like = "%F%"
np_like_1 = "%NP%"
np_like_2 = "%Np%"

# Optional final shapefile export
export_final_shp = True
out_export_folder = r"D:\Downloaded Shapes\NHD\Washington\exports"
out_export_name = "WABuffers_merged.shp"

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
def log(msg):
    print(msg)
    try:
        arcpy.AddMessage(msg)
    except Exception:
        pass


def ensure_field_exists(fc_path, field_name):
    fields = [f.name for f in arcpy.ListFields(fc_path)]
    if field_name not in fields:
        raise RuntimeError(f"Missing field '{field_name}' in {fc_path}")


def segment_half_away_from_confluence(geom, start_measure):
    length = geom.length
    if length <= 0:
        return None

    half_len = 0.5 * length
    d0 = start_measure
    d1 = length - start_measure
    desired_len = min(half_len, max(d0, d1))

    if d0 > d1:
        start = max(0.0, start_measure - desired_len)
        end = start_measure
    else:
        start = start_measure
        end = min(length, start_measure + desired_len)

    if end <= start:
        center = length / 2.0
        start = max(0.0, center - half_len / 2.0)
        end = min(length, center + half_len / 2.0)

    return geom.segmentAlongLine(start, end, use_percentage=False)


def to_point_geometry(any_geom, sr):
    if any_geom is None:
        return None
    try:
        if any_geom.type.lower() == "point":
            p = any_geom.firstPoint
        else:
            p = any_geom.centroid
    except Exception:
        p = any_geom.firstPoint
    if p is None:
        return None
    return arcpy.PointGeometry(arcpy.Point(p.X, p.Y), sr)


def is_endpoint(line_geom, pt_geom, eps=1e-6):
    if line_geom is None or pt_geom is None:
        return True

    sr = line_geom.spatialReference or arcpy.SpatialReference(target_sr_wkid)
    pt = to_point_geometry(pt_geom, sr)
    if pt is None:
        return True

    qpt, d_along, frac, lat_off = line_geom.queryPointAndDistance(pt, use_percentage=False)
    L = line_geom.length
    if L is None or L <= 0:
        return True
    return (d_along <= eps) or (abs(L - d_along) <= eps)


def field_delim_for(workspace, field):
    try:
        return arcpy.AddFieldDelimiters(workspace, field)
    except Exception:
        return field


def try_remove_mz(fc):
    """Attempt RemoveM/RemoveZ if available; skip gracefully otherwise."""
    for tool_name in ("RemoveM", "RemoveZ"):
        func = getattr(arcpy.management, tool_name, None)
        if func:
            try:
                func(fc)
                log(f"{tool_name} completed on {fc}.")
            except Exception as e:
                log(f"{tool_name} warning on {fc}: {e}")
        else:
            log(f"{tool_name} tool not available in this ArcGIS Pro environment; skipping.")


# ---------------------------------------------------------------------
# Workspace & CRS prep
# ---------------------------------------------------------------------
if not arcpy.Exists(out_gdb):
    folder = os.path.dirname(out_gdb)
    gdb_name = os.path.basename(out_gdb)
    arcpy.management.CreateFileGDB(folder, gdb_name)

scratch_gdb = arcpy.env.scratchGDB if arcpy.env.scratchGDB else out_gdb
target_sr = arcpy.SpatialReference(target_sr_wkid)

# ---------------------------------------------------------------------
# PRE-FLIGHT IMPORT to GDB (explicit name -> avoids Describe/path problems)
# ---------------------------------------------------------------------
log(f"Preflight import starting for: {flowlines_shp}")

# Repair geometry on the shapefile
try:
    arcpy.management.RepairGeometry(flowlines_shp, delete_null="DELETE_NULL")
    log("RepairGeometry complete on source shapefile.")
except Exception as e:
    log(f"RepairGeometry warning: {e}")

# Add spatial index to the shapefile
try:
    arcpy.management.AddSpatialIndex(flowlines_shp)
    log("AddSpatialIndex complete on source shapefile.")
except Exception as e:
    log(f"AddSpatialIndex warning: {e}")

# Explicit import name inside GDB
import_fc = os.path.join(out_gdb, "flowlines_src")
if not arcpy.Exists(import_fc):
    log(f"Importing shapefile to GDB as: {import_fc}")
    t0 = time.time()
    arcpy.conversion.FeatureClassToFeatureClass(
        in_features=flowlines_shp,
        out_path=out_gdb,
        out_name="flowlines_src",
    )
    log(f"FeatureClassToFeatureClass (import) took {time.time() - t0:.2f}s")
else:
    log("Import feature class already exists; reusing flowlines_src.")

# Optional: log hasM/hasZ after import
desc = arcpy.Describe(import_fc)
has_m = getattr(desc, "hasM", False)
has_z = getattr(desc, "hasZ", False)
log(f"Imported feature class hasM={has_m}, hasZ={has_z}")

# Optional: remove M/Z if tools are available
try_remove_mz(import_fc)

# Project the *GDB* feature class
flowlines_proj = os.path.join(scratch_gdb, "flowlines_projected")
log(f"Projecting to WKID {target_sr_wkid}: {flowlines_proj}")
t0 = time.time()
arcpy.management.Project(import_fc, flowlines_proj, target_sr)
log(f"Project took {time.time() - t0:.2f}s")
arcpy.env.outputCoordinateSystem = target_sr

# ---------------------------------------------------------------------
# Build layers for S, F, and Np (LIKE filters for robustness)
# ---------------------------------------------------------------------
workspace_for_sql = os.path.dirname(flowlines_proj)
type_f_delim = field_delim_for(workspace_for_sql, type_field)
s_lyr = "lyr_S"
f_lyr = "lyr_F"
np_lyr = "lyr_NP"
arcpy.management.MakeFeatureLayer(flowlines_proj, s_lyr, f"{type_f_delim} LIKE '{s_like}'")
arcpy.management.MakeFeatureLayer(flowlines_proj, f_lyr, f"{type_f_delim} LIKE '{f_like}'")
arcpy.management.MakeFeatureLayer(
    flowlines_proj,
    np_lyr,
    f"{type_f_delim} LIKE '{np_like_1}' OR {type_f_delim} LIKE '{np_like_2}'",
)

# ---------------------------------------------------------------------
# 1) Type S & F buffers = (half bankfull + 50 ft)
# ---------------------------------------------------------------------
sf_merge = os.path.join(scratch_gdb, "SF_lines_merge")
arcpy.management.Merge([s_lyr, f_lyr], sf_merge)
sf_lines_single_fc = os.path.join(scratch_gdb, "SF_lines_singlepart")
arcpy.management.MultipartToSinglepart(sf_merge, sf_lines_single_fc)
ensure_field_exists(sf_lines_single_fc, bankfull_field)
sf_buffers_fc = os.path.join(out_gdb, "buf_SF_halfBankfull_plus_50ft")
buf_field = "BUF_FT"
if buf_field not in [f.name for f in arcpy.ListFields(sf_lines_single_fc)]:
    arcpy.management.AddField(sf_lines_single_fc, buf_field, "DOUBLE")
with arcpy.da.UpdateCursor(sf_lines_single_fc, [bankfull_field, buf_field]) as ucur:
    for bfw, _ in ucur:
        bfw_val = float(bfw) if bfw is not None else 0.0
        half_width = bfw_val if is_half_width else (bfw_val * 0.5)
        ucur.updateRow([bfw_val, half_width + 50.0])
arcpy.analysis.Buffer(sf_lines_single_fc, sf_buffers_fc, buf_field, method="PLANAR")
sf_diss_fc = os.path.join(scratch_gdb, "SF_dissolved")
arcpy.management.Dissolve(sf_lines_single_fc, sf_diss_fc)

# ---------------------------------------------------------------------
# 2) Type Np: ONLY connected to S/F; first half AWAY from confluence; buffer (half+50)
# ---------------------------------------------------------------------
np_lines_copy_fc = os.path.join(scratch_gdb, "NP_lines_copy")
arcpy.management.CopyFeatures(np_lyr, np_lines_copy_fc)
if "ORIG_ID_NP" not in [f.name for f in arcpy.ListFields(np_lines_copy_fc)]:
    arcpy.management.AddField(np_lines_copy_fc, "ORIG_ID_NP", "LONG")
    with arcpy.da.UpdateCursor(np_lines_copy_fc, ["ORIG_ID_NP", "OID@"]) as ucur:
        for _, oid0 in ucur:
            ucur.updateRow([oid0, oid0])

np_lines_single_fc = os.path.join(scratch_gdb, "NP_lines_singlepart")
arcpy.management.MultipartToSinglepart(np_lines_copy_fc, np_lines_single_fc)
ensure_field_exists(np_lines_single_fc, bankfull_field)
if "ORIG_ID_NP" not in [f.name for f in arcpy.ListFields(np_lines_single_fc)]:
    np_lines_single_join_fc = os.path.join(scratch_gdb, "NP_lines_singlepart_join")
    arcpy.analysis.SpatialJoin(
        target_features=np_lines_single_fc,
        join_features=np_lines_copy_fc,
        out_feature_class=np_lines_single_join_fc,
        join_operation="JOIN_ONE_TO_ONE",
        join_type="KEEP_ALL",
        match_option="INTERSECT",
    )
    np_lines_single_fc = np_lines_single_join_fc
    ensure_field_exists(np_lines_single_fc, "ORIG_ID_NP")

near_table = os.path.join(scratch_gdb, "NP_to_SF_near")
arcpy.analysis.GenerateNearTable(
    in_features=np_lines_single_fc,
    near_features=sf_diss_fc,
    out_table=near_table,
    search_radius=f"{snap_tol_ft} Feet",
    location="LOCATION",
    angle="NO_ANGLE",
    closest="ALL",
    method="PLANAR",
)

near_index = {}
with arcpy.da.SearchCursor(near_table, ["IN_FID", "NEAR_FID", "NEAR_DIST", "NEAR_X", "NEAR_Y"]) as sc:
    for in_fid, near_fid, near_dist, nx, ny in sc:
        if near_dist is None:
            continue
        if float(near_dist) <= snap_tol_ft:
            if in_fid not in near_index or float(near_dist) < near_index[in_fid][2]:
                near_index[in_fid] = (float(nx), float(ny), float(near_dist))

np_half_fc = os.path.join(out_gdb, "NP_first_half_segments")
arcpy.management.CreateFeatureclass(
    out_gdb,
    os.path.basename(np_half_fc),
    "POLYLINE",
    spatial_reference=target_sr,
)
if bankfull_field not in [f.name for f in arcpy.ListFields(np_half_fc)]:
    arcpy.management.AddField(np_half_fc, bankfull_field, "DOUBLE")

fields_search = ["OID@", "SHAPE@", bankfull_field]
fields_insert = ["SHAPE@", bankfull_field]
skipped_np = 0
with arcpy.da.InsertCursor(np_half_fc, fields_insert) as icur, arcpy.da.SearchCursor(
    np_lines_single_fc, fields_search
) as scur:
    for oid, geom, bfw in scur:
        if oid not in near_index:
            skipped_np += 1
            continue
        nx, ny, _ = near_index[oid]
        bfw_val = float(bfw) if bfw is not None else 0.0
        near_pt = arcpy.PointGeometry(arcpy.Point(nx, ny), target_sr)
        qpt, distAlong, frac, latOff = geom.queryPointAndDistance(near_pt, use_percentage=False)
        seg = segment_half_away_from_confluence(geom, float(distAlong))
        if seg is None or seg.length <= 0:
            skipped_np += 1
            continue
        icur.insertRow([seg, bfw_val])
log(f"NP first-half segments created: {int(arcpy.management.GetCount(np_half_fc)[0])}, skipped: {skipped_np}")

np_half_buffers_fc = os.path.join(out_gdb, "buf_NP_first_half_halfBankfull_plus_50ft")
buf_field_np = "BUF_FT"
if buf_field_np not in [f.name for f in arcpy.ListFields(np_half_fc)]:
    arcpy.management.AddField(np_half_fc, buf_field_np, "DOUBLE")
with arcpy.da.UpdateCursor(np_half_fc, [bankfull_field, buf_field_np]) as ucur:
    for bfw_val, _ in ucur:
        bfw_val = float(bfw_val) if bfw_val is not None else 0.0
        half_width = bfw_val if is_half_width else (bfw_val * 0.5)
        ucur.updateRow([bfw_val, half_width + 50.0])
arcpy.analysis.Buffer(np_half_fc, np_half_buffers_fc, buf_field_np, method="PLANAR")

# ---------------------------------------------------------------------
# 3) Type Np intersections (Western WA): true intersections only
#     (Improved: FeatureVerticesToPoints INTERSECTION)
# ---------------------------------------------------------------------
np_x_pts_raw_fc = os.path.join(scratch_gdb, "NP_intersections_pts_raw")
arcpy.management.FeatureVerticesToPoints(np_lines_single_fc, np_x_pts_raw_fc, "INTERSECTION")

# Deduplicate identical geometry points
arcpy.management.DeleteIdentical(np_x_pts_raw_fc, ["Shape"])

# Layer to select touching Np lines per point
np_lines_layer = "np_lines_layer"
arcpy.management.MakeFeatureLayer(np_lines_single_fc, np_lines_layer)

# Output: only true intersections
np_x_pts_true_fc = os.path.join(scratch_gdb, "NP_intersections_pts_true")
arcpy.management.CreateFeatureclass(
    scratch_gdb,
    os.path.basename(np_x_pts_true_fc),
    "POINT",
    spatial_reference=target_sr,
)
if "PT_ID" not in [f.name for f in arcpy.ListFields(np_x_pts_true_fc)]:
    arcpy.management.AddField(np_x_pts_true_fc, "PT_ID", "LONG")

kept = 0
checked = 0
with arcpy.da.InsertCursor(np_x_pts_true_fc, ["PT_ID", "SHAPE@"]) as ic_true, arcpy.da.SearchCursor(
    np_x_pts_raw_fc, ["OID@", "SHAPE@"]
) as sc_pts:
    for pt_oid, pt_geom in sc_pts:
        checked += 1
        arcpy.management.SelectLayerByLocation(
            np_lines_layer,
            "INTERSECT",
            pt_geom,
            selection_type="NEW_SELECTION",
        )
        distinct_line_oids = set()
        distinct_orig_ids = set()
        interior_count = 0
        with arcpy.da.SearchCursor(np_lines_layer, ["OID@", "ORIG_ID_NP", "SHAPE@"]) as sc_lines:
            for line_oid, orig_id, line_geom in sc_lines:
                distinct_line_oids.add(line_oid)
                if orig_id is not None:
                    distinct_orig_ids.add(int(orig_id))
                if not is_endpoint(line_geom, pt_geom):
                    interior_count += 1
        if (len(distinct_line_oids) >= 2) and (len(distinct_orig_ids) >= 2) and (interior_count >= 1):
            kept += 1
            ic_true.insertRow([kept, to_point_geometry(pt_geom, target_sr)])
log(f"Checked Np intersection candidates: {checked}, kept true intersections: {kept}")
log(f"True intersection points count (to buffer): {int(arcpy.management.GetCount(np_x_pts_true_fc)[0])}")

# Buffer true intersections to fixed Western WA radius (56 ft)
np_x_buffers_fc = os.path.join(out_gdb, "buf_NP_intersections_56ft")
arcpy.analysis.Buffer(np_x_pts_true_fc, np_x_buffers_fc, f"{intersection_radius_ft} Feet", method="PLANAR")

# ---------------------------------------------------------------------
# 4) Optional merged layer (all buffers)
# ---------------------------------------------------------------------
merged_fc = os.path.join(out_gdb, "buf_all_merged")
arcpy.management.Merge([sf_buffers_fc, np_half_buffers_fc, np_x_buffers_fc], merged_fc)
log("Completed.")
log(f" S/F buffers (half bankfull + 50 ft): {sf_buffers_fc}")
log(f" Np first-half segments: {np_half_fc}")
log(f" Np first-half buffers (half+50): {np_half_buffers_fc}")
log(f" Np intersections (true points): {np_x_pts_true_fc}")
log(f" Np intersections buffers (56'): {np_x_buffers_fc}")
log(f" Merged buffers: {merged_fc}")
log(f" Merged count: {int(arcpy.management.GetCount(merged_fc)[0])}")

# ---------------------------------------------------------------------
# 5) Safe export to shapefile (optional)
# ---------------------------------------------------------------------
if export_final_shp:
    try:
        arcpy.management.RepairGeometry(merged_fc, delete_null="DELETE_NULL")
        try_remove_mz(merged_fc)
        if not os.path.isdir(out_export_folder):
            os.makedirs(out_export_folder)
        out_shp = os.path.join(out_export_folder, out_export_name)
        log(f"Exporting final shapefile to: {out_shp}")
        keep_fields = ["BUF_FT", bankfull_field]
        fm = arcpy.FieldMappings()
        fm.addTable(merged_fc)
        for fld in [f.name for f in arcpy.ListFields(merged_fc) if f.name not in keep_fields + ["Shape"]]:
            idx = fm.findFieldMapIndex(fld)
            if idx != -1:
                fm.removeFieldMap(idx)

        def rename_field(old_name, new_name):
            idx = fm.findFieldMapIndex(old_name)
            if idx != -1:
                fmap = fm.getFieldMap(idx)
                ofld = fmap.outputField
                ofld.name = new_name  # <=10 chars for .shp
                fmap.outputField = ofld
                fm.replaceFieldMap(idx, fmap)

        rename_field(bankfull_field, "BANKFULL")
        rename_field("BUF_FT", "BUF_FT")
        t0 = time.time()
        arcpy.conversion.FeatureClassToFeatureClass(
            in_features=merged_fc,
            out_path=out_export_folder,
            out_name=out_export_name,
            field_mapping=fm,
        )
        log(f"FeatureClassToFeatureClass export completed in {time.time() - t0:.2f}s")
        try:
            arcpy.management.AddSpatialIndex(out_shp)
        except Exception:
            pass
        log("Final shapefile written successfully.")
        log(out_shp)
    except Exception as e:
        log(f"Export to shapefile failed: {e}")

# End-of-script status
log("Workflow completed.")
if export_final_shp:
    log("Final shapefile export attempted.")
