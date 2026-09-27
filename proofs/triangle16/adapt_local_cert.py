"""STRUCTURAL adapter (2026-09-26, lead): LOCAL.md did not fix the JSON layout, so the blind checker (check_local_blind.py,
written from LOCAL.md alone) expects cert['clusters'] = [{id, frame, rattler, pairs, walls:[{point,kind}], lambda,
y:[{coord,sign,values}], ...}] while local_cert.py writes parallel lists 'summary'/'data'. This script ONLY renames and
re-nests fields; it never changes a value. Mappings: id = summary[i]['cluster']; frame = data['frame_uv']; rattler = the single
entry of data['rattlers_uv']; walls [a, 'u'|'v'|'w'] -> {point: a, kind: 'u'|'v'|'s'} (the checker names the slanted side u+v=1
's'); y[2i + (0 if sign=+1 else 1)] -> {coord: i, sign: +-1, values} (local_cert.py loops i in range(30), sign in (1, -1)).
Every other data field is passed through unchanged. A wrong mapping can only make the checker FAIL (it recomputes every
identity from the stored rationals), never PASS.
usage: python adapt_local_cert.py out/local_cert16.json out/local_cert16_for_blind.json"""
import json, sys
src, dst = sys.argv[1], sys.argv[2]
C = json.load(open(src))
assert len(C["summary"]) == len(C["data"]) == 3
out = []
for s, d in zip(C["summary"], C["data"]):
    assert len(d["rattlers_uv"]) == 1 and len(d["y"]) == 60 and len(d["frame_uv"]) == 15
    item = {k: v for k, v in d.items() if k not in ("frame_uv", "rattlers_uv", "walls", "y")}
    item["id"] = s["cluster"]
    item["frame"] = d["frame_uv"]
    (rk, rv), = d["rattlers_uv"].items(); item["rattler"] = rv; item["rattler_label"] = int(rk)
    item["walls"] = [{"point": a, "kind": {"u": "u", "v": "v", "w": "s"}[k]} for a, k in d["walls"]]
    item["y"] = [{"coord": n // 2, "sign": 1 if n % 2 == 0 else -1, "values": v} for n, v in enumerate(d["y"])]
    out.append(item)
json.dump({"clusters": out}, open(dst, "w"), indent=0)
if len(sys.argv) > 4:                     # boxes: frame = the 15 masked points (label order), as (A, B) of the barycentric (A, B, C)
    Bx = json.load(open(sys.argv[3]))
    for c in Bx["clusters"]:
        assert len(c["points"]) == len(c["frame"]) == 16 and sum(c["frame"]) == 15
        c["frame_mask"] = c["frame"]; c["frame"] = [[pt[0], pt[1]] for pt, keep in zip(c["points"], c["frame_mask"]) if keep]
    json.dump(Bx, open(sys.argv[4], "w"), indent=0); print("adapted boxes ->", sys.argv[4])
print("adapted", len(out), "clusters ->", dst)
