#!/usr/bin/env python3
import argparse,hashlib,importlib.util,json,math,pathlib,sys
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--package-dir",required=True);a=ap.parse_args();root=pathlib.Path(a.package_dir);man=json.loads((root/"MANIFEST.json").read_text())
 for x in man["files"]:
  p=root/x["path"]
  if sha(p)!=x["sha256"] or p.stat().st_size!=x["bytes"]:raise ValueError("manifest mismatch: "+x["path"])
 s=importlib.util.spec_from_file_location("portable",root/"portable_hgb_predictor.py");m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
 f=json.loads((root/"SYNTHETIC_EQUIVALENCE_FIXTURES.json").read_text());mp=root/(f["model_id"]+".portable.json");doc=m.load_model(mp,sha(mp));rows=[[float.fromhex(v) for v in r] for r in f["rows_float_hex"]];obs=m.predict_positive_proba(doc,rows,f["feature_order"]);exp=[float.fromhex(v) for v in f["expected_positive_probability_float_hex"]];t=float.fromhex(doc["decision_threshold"])
 for i,(x,y) in enumerate(zip(obs,exp)):
  if not math.isclose(x,y,rel_tol=1e-12,abs_tol=1e-12) or int(x>=t)!=f["expected_calls"][i]:raise ValueError("self-test mismatch row "+str(i))
 print("PASS: portable documented-EdU G1 synthetic self-test")
if __name__=="__main__":
 try:main()
 except Exception as e:print("STOP: "+str(e),file=sys.stderr);sys.exit(2)
