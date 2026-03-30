diff --git a/WABuffer.py b/WABuffer.py
index bf83228dcdfa2428e4cfedf98a3b53075370f911..94544920c19cea771dd6d7a916a01f99acbbb800 100644
--- a/WABuffer.py
+++ b/WABuffer.py
@@ -1,58 +1,104 @@
 # SMZWWA_fixed.py
 # ArcGIS Pro / ArcPy 3.x (Python 3.9+)
 # Western WA workflow:
 # 1) Type S & F buffers: (half bankfull + 50 ft), centerline-based
 # 2) Type Np: ONLY Np connected to S/F; take FIRST HALF away from confluence toward upstream; buffer (half bankfull + 50 ft)
 # 3) Type Np intersections (Western WA): fixed 56-ft buffers (true intersections only: ≥2 distinct Np streams and ≥1 interior contact)
 # 4) Optional export of final merged buffers to shapefile (safe field mapping, optional RemoveZ/RemoveM, spatial index)
 import arcpy
+import argparse
 import os
+import re
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
 target_sr_wkid = 2913          # WA North (ftUS); use 2914 for WA South
 snap_tol_ft = 10.0
 intersection_radius_ft = 56.0  # fixed Western WA site radius
 is_half_width = False
 out_gdb = r"D:\Downloaded Shapes\NHD\Washington\WABuffers.gdb"
-# LIKE patterns
-s_like = "%S%"
-f_like = "%F%"
-np_like_1 = "%NP%"
-np_like_2 = "%Np%"
+# Canonical type labels used by this workflow
+TYPE_S = "S"
+TYPE_F = "F"
+TYPE_NP = "NP"
 # Optional final shapefile export
 export_final_shp = True
 out_export_folder = r"D:\Downloaded Shapes\NHD\Washington\exports"
 out_export_name = "WABuffers_merged.shp"
+def parse_args():
+    """
+    Optional CLI overrides for running directly in VS Code terminal.
+    Any values not supplied on CLI keep the defaults above.
+    """
+    parser = argparse.ArgumentParser(
+        description="Build Western WA SMZ/no-harvest buffers from hydro flowlines."
+    )
+    parser.add_argument("--flowlines", dest="flowlines_shp", help="Input flowlines shapefile/feature class path.")
+    parser.add_argument("--type-field", dest="type_field", help="Field that stores stream type code (e.g., FP_REF_ID).")
+    parser.add_argument("--bankfull-field", dest="bankfull_field", help="Field containing full bankfull width in feet.")
+    parser.add_argument("--out-gdb", dest="out_gdb", help="Output file geodatabase path.")
+    parser.add_argument("--target-wkid", dest="target_sr_wkid", type=int, help="Projected output WKID (2913 WA North / 2914 WA South).")
+    parser.add_argument("--snap-tol-ft", dest="snap_tol_ft", type=float, help="Near-table snap tolerance in feet.")
+    parser.add_argument("--intersection-radius-ft", dest="intersection_radius_ft", type=float, help="NP intersection fixed buffer radius in feet.")
+    parser.add_argument("--is-half-width", dest="is_half_width", action="store_true", help="Set when bankfull field is already half-width.")
+    parser.add_argument("--no-export-shp", dest="no_export_shp", action="store_true", help="Skip final shapefile export.")
+    parser.add_argument("--export-folder", dest="out_export_folder", help="Folder for final exported shapefile.")
+    parser.add_argument("--export-name", dest="out_export_name", help="Filename for final exported shapefile.")
+    args, _unknown = parser.parse_known_args()
+    return args
+
+cli = parse_args()
+if cli.flowlines_shp:
+    flowlines_shp = cli.flowlines_shp
+if cli.type_field:
+    type_field = cli.type_field
+if cli.bankfull_field:
+    bankfull_field = cli.bankfull_field
+if cli.out_gdb:
+    out_gdb = cli.out_gdb
+if cli.target_sr_wkid:
+    target_sr_wkid = cli.target_sr_wkid
+if cli.snap_tol_ft is not None:
+    snap_tol_ft = cli.snap_tol_ft
+if cli.intersection_radius_ft is not None:
+    intersection_radius_ft = cli.intersection_radius_ft
+if cli.is_half_width:
+    is_half_width = True
+if cli.no_export_shp:
+    export_final_shp = False
+if cli.out_export_folder:
+    out_export_folder = cli.out_export_folder
+if cli.out_export_name:
+    out_export_name = cli.out_export_name
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
@@ -71,125 +117,165 @@ def to_point_geometry(any_geom, sr):
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
+def canonical_stream_type(raw_value):
+    """
+    Normalize stream type text to one of: S, F, NP, or None.
+    Handles mixed case and common separators (space, -, _, /).
+    """
+    if raw_value is None:
+        return None
+    txt = str(raw_value).strip().upper()
+    if not txt:
+        return None
+    # Normalize separators to spaces for token checks.
+    normalized = re.sub(r"[\-_/]+", " ", txt)
+    tokens = [t for t in re.split(r"\s+", normalized) if t]
+    if any(t.startswith("NP") for t in tokens):
+        return TYPE_NP
+    if any(t.startswith("S") for t in tokens):
+        return TYPE_S
+    if any(t.startswith("F") for t in tokens):
+        return TYPE_F
+    return None
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
+if not arcpy.Exists(flowlines_shp):
+    raise RuntimeError(
+        "Input flowlines path does not exist. "
+        "Set 'flowlines_shp' in the script or pass --flowlines in VS Code terminal."
+    )
+log("Running WABuffer with configuration:")
+log(f" flowlines_shp={flowlines_shp}")
+log(f" out_gdb={out_gdb}")
+log(f" type_field={type_field}, bankfull_field={bankfull_field}")
+log(f" target_sr_wkid={target_sr_wkid}, snap_tol_ft={snap_tol_ft}, intersection_radius_ft={intersection_radius_ft}")
+log(f" is_half_width={is_half_width}, export_final_shp={export_final_shp}")
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
         out_name="flowlines_src"
     )
     log(f"FeatureClassToFeatureClass (import) took {time.time()-t0:.2f}s")
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
 log(f"Project took {time.time()-t0:.2f}s")
 arcpy.env.outputCoordinateSystem = target_sr
+# Prepare a normalized type field so selection is robust/case-insensitive
+flowlines_typed = os.path.join(scratch_gdb, "flowlines_projected_typed")
+arcpy.management.CopyFeatures(flowlines_proj, flowlines_typed)
+wa_type_field = "WA_TYPE"
+if wa_type_field not in [f.name for f in arcpy.ListFields(flowlines_typed)]:
+    arcpy.management.AddField(flowlines_typed, wa_type_field, "TEXT", field_length=8)
+with arcpy.da.UpdateCursor(flowlines_typed, [type_field, wa_type_field]) as ucur:
+    for type_raw, _ in ucur:
+        ucur.updateRow([type_raw, canonical_stream_type(type_raw)])
 # ---------------------------------------------------------------------
-# Build layers for S, F, and Np (LIKE filters for robustness)
+# Build layers for S, F, and Np from canonicalized type values
 # ---------------------------------------------------------------------
-workspace_for_sql = os.path.dirname(flowlines_proj)
-type_f_delim = field_delim_for(workspace_for_sql, type_field)
+workspace_for_sql = os.path.dirname(flowlines_typed)
+wa_type_delim = field_delim_for(workspace_for_sql, wa_type_field)
 s_lyr = "lyr_S"
 f_lyr = "lyr_F"
 np_lyr = "lyr_NP"
-arcpy.management.MakeFeatureLayer(flowlines_proj, s_lyr, f"{type_f_delim} LIKE '{s_like}'")
-arcpy.management.MakeFeatureLayer(flowlines_proj, f_lyr, f"{type_f_delim} LIKE '{f_like}'")
-arcpy.management.MakeFeatureLayer(flowlines_proj, np_lyr, f"{type_f_delim} LIKE '{np_like_1}' OR {type_f_delim} LIKE '{np_like_2}'")
+arcpy.management.MakeFeatureLayer(flowlines_typed, s_lyr, f"{wa_type_delim} = '{TYPE_S}'")
+arcpy.management.MakeFeatureLayer(flowlines_typed, f_lyr, f"{wa_type_delim} = '{TYPE_F}'")
+arcpy.management.MakeFeatureLayer(flowlines_typed, np_lyr, f"{wa_type_delim} = '{TYPE_NP}'")
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
