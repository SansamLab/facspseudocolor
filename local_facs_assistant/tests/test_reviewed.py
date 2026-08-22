"""Milestone 3 tests using only clearly labeled SYNTHETIC fixtures and mocks."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import os

from jsonschema import Draft202012Validator, ValidationError

from local_facs_assistant.facs_tools import IntakeError
from local_facs_assistant.reviewed import (
    CONFIRMATION_STATEMENT, METHOD_AUTHORIZATION, QMD_NAME, REVIEWED_NAME,
    _control_relationships, _hash_regular_file, _read_bounded_once,
    confirm_proposal, generate_scaffold, save_reviewed,
    validate_identity,
)


class ReviewedConfigTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="SYNTHETIC_REVIEWED_")
        self.root = Path(self.temp.name)
        for name in ("SYNTHETIC.fcs", "SYNTHETIC_sample.fcs", "SYNTHETIC.wsp"):
            (self.root / name).write_text("SYNTHETIC", encoding="utf-8")
        self.proposal = {
            "schema_version": "2.0",
            "experiment": {"directory": self.root.name, "title": None, "biological_replicates": None},
            "inputs": {"fcs_files": ["SYNTHETIC.fcs","SYNTHETIC_sample.fcs"], "workspace": "SYNTHETIC.wsp"},
            "sample_mapping": [
                {"file":"SYNTHETIC.fcs","condition":"SYNTHETIC background","time":"not applicable","role":"matched_background_control","biological_replicate":1,"provenance":"user_supplied","confirmed":False,"evidence":["SYNTHETIC.fcs"]},
                {"file":"SYNTHETIC_sample.fcs","condition":"SYNTHETIC condition","time":"0 h","role":"experimental_sample","biological_replicate":1,"provenance":"user_supplied","confirmed":False,"evidence":["SYNTHETIC_sample.fcs"]}],
            "channels": [
                {"category":"DNA","feature":"SYNTHETIC DNA","detector":"FL2-A","label":"PE-A","provenance":"user_supplied","confirmed":False,"evidence":["SYNTHETIC.fcs","SYNTHETIC_sample.fcs"]},
                {"category":"POI","feature":"SYNTHETIC POI","detector":"FL4-A","label":"APC-A","provenance":"user_supplied","confirmed":False,"evidence":["SYNTHETIC.fcs","SYNTHETIC_sample.fcs"]},
            ],
            "analyses": [{"name":"SYNTHETIC POI vs DNA","analysis_type":"poi_vs_dna","target_feature":"SYNTHETIC POI","dna_feature":"SYNTHETIC DNA","population":"Single Cells","provenance":"user_supplied","confirmed":False}],
            "recorded_details": [{"workspace":"SYNTHETIC.wsp","layout":"SYNTHETIC layout","text":"SYNTHETIC note"}],
            "model_review": {"provenance":"model_advisory","status":"no_additional_uncertainty","flags":[]},
            "authorization": {"sample_mapping_confirmed":False,"channel_mapping_confirmed":False,"analysis_selection_confirmed":False,"analysis_authorized":False},
            "uncertainties": [], "evidence": ["SYNTHETIC.fcs","SYNTHETIC_sample.fcs", "SYNTHETIC.wsp"],
        }
        self.proposal_path = self.root / "SYNTHETIC_proposal.json"
        self.proposal_path.write_text(json.dumps(self.proposal), encoding="utf-8")
        self.ledger = [
            {"tool":"inventory_experiment","arguments":{},"output":{"experiment_root":self.root.name,"files":[
                {"path":"SYNTHETIC.fcs","suffix":".fcs","bytes":9},{"path":"SYNTHETIC_sample.fcs","suffix":".fcs","bytes":9}, {"path":"SYNTHETIC.wsp","suffix":".wsp","bytes":9}]}},
            {"tool":"inspect_fcs_metadata","arguments":{"path":"SYNTHETIC.fcs"},"output":{"file":"SYNTHETIC.fcs","event_count":1,"channels":[
                {"index":1,"detector":"FL2-A","stain":"PE-A"},{"index":2,"detector":"FL4-A","stain":"APC-A"}],"metadata":{}}},
            {"tool":"inspect_fcs_metadata","arguments":{"path":"SYNTHETIC_sample.fcs"},"output":{"file":"SYNTHETIC_sample.fcs","event_count":1,"channels":[
                {"index":1,"detector":"FL2-A","stain":"PE-A"},{"index":2,"detector":"FL4-A","stain":"APC-A"}],"metadata":{}}},
            {"tool":"inspect_wsp","arguments":{"path":"SYNTHETIC.wsp"},"output":{"workspace":"SYNTHETIC.wsp","samples":[
                {"file":"SYNTHETIC.fcs","workspace_sample_id":"1","gate_names":["Single Cells"],"gates":[{"name":"Single Cells","path":["Single Cells"]}]},
                {"file":"SYNTHETIC_sample.fcs","workspace_sample_id":"2","gate_names":["Single Cells"],"gates":[{"name":"Single Cells","path":["Single Cells"]}]}]}},
            {"tool":"extract_layout_text","arguments":{"path":"SYNTHETIC.wsp"},"output":{"workspace":"SYNTHETIC.wsp","layouts":[{"layout":"SYNTHETIC layout","text":["SYNTHETIC note"]}]}},
        ]

    def tearDown(self): self.temp.cleanup()

    def _challenge(self): return hashlib.sha256(self.proposal_path.read_bytes()).hexdigest()
    def _relationship(self, **changes):
        value={"control_file":"SYNTHETIC.fcs","relationship":"matched_background_control",
            "applies_to_samples":["SYNTHETIC_sample.fcs"],"applies_to_analyses":[self.proposal["analyses"][0]["name"]],
            "applies_to_features":[self.proposal["analyses"][0]["target_feature"]]}
        value.update(changes); return json.dumps(value)

    def _confirm(self):
        with patch("local_facs_assistant.reviewed.preflight_inspect", return_value=(self.ledger, "SYNTHETIC")):
            return confirm_proposal(self.root,self.proposal_path,CONFIRMATION_STATEMENT,self._challenge(),[self._relationship()])

    def test_confirmation_is_explicit_complete_and_methods_stay_false(self):
        reviewed = self._confirm()
        self.assertTrue(all(row["confirmed"] is True for row in reviewed["sample_mapping"] + reviewed["channels"] + reviewed["analyses"]))
        self.assertTrue(all(value is True for key, value in reviewed["confirmation"].items() if key.endswith("_confirmed")))
        self.assertEqual(reviewed["method_authorization"], METHOD_AUTHORIZATION)
        schema=json.loads((Path(__file__).parents[1]/"schemas/reviewed-config.schema.json").read_text())
        Draft202012Validator(schema).validate(reviewed)
        with self.assertRaises(IntakeError):
            confirm_proposal(self.root,self.proposal_path,"I confirm something else",self._challenge(),[self._relationship()])
        changed=copy.deepcopy(reviewed); changed["confirmation"]["analysis_selection_confirmed"]=False
        with self.assertRaises(ValidationError): Draft202012Validator(schema).validate(changed)

    def test_confirmation_preserves_proposal_values_without_layout_enrichment(self):
        self.proposal["channels"][1]["feature"]="SYNTHETIC phospho-FOXM1"
        self.proposal["analyses"][0]["name"]="SYNTHETIC pFOXM1 vs DNA"
        self.proposal["analyses"][0]["target_feature"]="SYNTHETIC phospho-FOXM1"
        self.proposal["recorded_details"][0]["text"]="SYNTHETIC note says phospho-FOXM1 T600"
        self.ledger[4]["output"]["layouts"][0]["text"]=["SYNTHETIC note says phospho-FOXM1 T600"]
        self.proposal_path.write_text(json.dumps(self.proposal),encoding="utf-8")
        reviewed=self._confirm()
        for field in ("experiment","inputs","recorded_details","uncertainties","evidence"):
            self.assertEqual(reviewed[field],self.proposal[field])
        for field in ("sample_mapping","channels","analyses"):
            expected=copy.deepcopy(self.proposal[field])
            for row in expected: row["confirmed"]=True
            self.assertEqual(reviewed[field],expected)
        self.assertEqual(reviewed["channels"][1]["feature"],"SYNTHETIC phospho-FOXM1")
        self.assertNotIn("T600",reviewed["channels"][1]["feature"])

    def test_confirmation_rejects_incomplete_claims(self):
        relationship=self._relationship()
        self.proposal["analyses"]=[]
        self.proposal_path.write_text(json.dumps(self.proposal),encoding="utf-8")
        with patch("local_facs_assistant.reviewed.preflight_inspect", return_value=(self.ledger,"SYNTHETIC")):
            with self.assertRaises(IntakeError): confirm_proposal(self.root,self.proposal_path,CONFIRMATION_STATEMENT,self._challenge(),[relationship])

    def test_control_relationship_scope_is_exact(self):
        self.assertEqual(len(_control_relationships([self._relationship()],self.proposal)),1)
        invalid=(
            self._relationship(applies_to_samples=[]),
            self._relationship(applies_to_samples=["SYNTHETIC_sample.fcs","UNKNOWN.fcs"]),
            self._relationship(applies_to_analyses=["UNKNOWN analysis"]),
            self._relationship(applies_to_features=["UNKNOWN feature"]),
        )
        for claim in invalid:
            with self.assertRaises(IntakeError): _control_relationships([claim],self.proposal)
        with self.assertRaises(IntakeError): _control_relationships([self._relationship(),self._relationship()],self.proposal)

    def test_multiple_controls_partition_cartesian_scope_without_overlap(self):
        proposal=copy.deepcopy(self.proposal)
        for file,role in (("SYNTHETIC_B.fcs","experimental_sample"),("SYNTHETIC_C.fcs","experimental_sample"),("SYNTHETIC_ctrl2.fcs","matched_background_control")):
            proposal["sample_mapping"].append({"file":file,"condition":"SYNTHETIC","time":"1 h","role":role,"biological_replicate":2,"provenance":"user_supplied","confirmed":False,"evidence":[file]})
        base={"relationship":"matched_background_control","applies_to_analyses":["SYNTHETIC POI vs DNA"],"applies_to_features":["SYNTHETIC POI"]}
        first=json.dumps({**base,"control_file":"SYNTHETIC.fcs","applies_to_samples":["SYNTHETIC_sample.fcs"]})
        second=json.dumps({**base,"control_file":"SYNTHETIC_ctrl2.fcs","applies_to_samples":["SYNTHETIC_B.fcs","SYNTHETIC_C.fcs"]})
        self.assertEqual(len(_control_relationships([first,second],proposal)),2)
        first_overlap=json.dumps({**base,"control_file":"SYNTHETIC.fcs","applies_to_samples":["SYNTHETIC_sample.fcs","SYNTHETIC_B.fcs"]})
        overlap=json.dumps({**base,"control_file":"SYNTHETIC_ctrl2.fcs","applies_to_samples":["SYNTHETIC_B.fcs","SYNTHETIC_C.fcs"]})
        with self.assertRaises(IntakeError): _control_relationships([first_overlap,overlap],proposal)
        gap=json.dumps({**base,"control_file":"SYNTHETIC_ctrl2.fcs","applies_to_samples":["SYNTHETIC_B.fcs"]})
        with self.assertRaises(IntakeError): _control_relationships([self._relationship(),gap],proposal)
        background_target=json.dumps({**base,"control_file":"SYNTHETIC.fcs","applies_to_samples":["SYNTHETIC_ctrl2.fcs"]})
        with self.assertRaises(IntakeError): _control_relationships([background_target],proposal)

    def test_no_follow_fd_read_rejects_path_replacement(self):
        replacement=self.root/"SYNTHETIC_replacement.json"
        replacement.write_bytes(self.proposal_path.read_bytes())
        real_read=os.read; replaced=False
        def replace_after_read(fd,size):
            nonlocal replaced
            data=real_read(fd,size)
            if data and not replaced:
                replaced=True; os.replace(replacement,self.proposal_path)
            return data
        with patch("local_facs_assistant.reviewed.os.read",side_effect=replace_after_read):
            with self.assertRaises(IntakeError): _read_bounded_once(self.proposal_path,2*1024*1024)

    def test_no_follow_input_hash_rejects_path_replacement(self):
        source=self.root/"SYNTHETIC_sample.fcs"; replacement=self.root/"SYNTHETIC_replacement.fcs"
        replacement.write_bytes(source.read_bytes())
        real_read=os.read; replaced=False
        def replace_after_read(fd,size):
            nonlocal replaced
            data=real_read(fd,size)
            if data and not replaced:
                replaced=True; os.replace(replacement,source)
            return data
        with patch("local_facs_assistant.reviewed.os.read",side_effect=replace_after_read):
            with self.assertRaises(IntakeError): _hash_regular_file(source,1024*1024)

    def test_digest_challenge_and_change_during_inspection_fail(self):
        with self.assertRaises(IntakeError):
            confirm_proposal(self.root,self.proposal_path,CONFIRMATION_STATEMENT,"0"*64,[self._relationship()])
        before={"SYNTHETIC.fcs":{"size_bytes":1,"sha256":"1"*64}}
        after={"SYNTHETIC.fcs":{"size_bytes":2,"sha256":"2"*64}}
        with patch("local_facs_assistant.reviewed._hash_inputs",side_effect=[before,after]), patch("local_facs_assistant.reviewed.preflight_inspect",return_value=(self.ledger,"SYNTHETIC")):
            with self.assertRaises(IntakeError): confirm_proposal(self.root,self.proposal_path,CONFIRMATION_STATEMENT,self._challenge(),[self._relationship()])

    def test_tampered_claim_fails_reinspection_and_nothing_is_written(self):
        changed=copy.deepcopy(self.proposal); changed["channels"][1]["detector"]="FABRICATED"
        self.proposal_path.write_text(json.dumps(changed), encoding="utf-8")
        with patch("local_facs_assistant.reviewed.preflight_inspect", return_value=(self.ledger, "SYNTHETIC")):
            with self.assertRaises(IntakeError): save_reviewed(self.root,self.proposal_path,CONFIRMATION_STATEMENT,self._challenge(),[self._relationship()])
        self.assertFalse((self.root/REVIEWED_NAME).exists())

    def test_generate_is_deterministic_escaped_and_has_method_barrier(self):
        injected='SYNTHETIC ":\nexecute: true #'
        self.proposal["analyses"][0]["name"]=injected
        self.proposal_path.write_text(json.dumps(self.proposal),encoding="utf-8")
        reviewed=self._confirm()
        reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed), encoding="utf-8")
        first=generate_scaffold(self.root, reviewed_path)
        bytes_first={path.name:path.read_bytes() for path in first}
        qmd=(self.root/QMD_NAME).read_text(encoding="utf-8")
        self.assertIn("authorization-barrier",qmd)
        self.assertLess(qmd.index("authorization-barrier"),qmd.index("facspseudocolor"))
        self.assertNotIn("quarto render",qmd); self.assertNotIn("export_flowjo",qmd)
        self.assertIn("render_authorized: false",qmd)
        for path in first: path.unlink()
        second=generate_scaffold(self.root, reviewed_path)
        self.assertEqual(bytes_first,{path.name:path.read_bytes() for path in second})
        identity=json.loads(next(path for path in second if path.suffix==".json" and path.name!=REVIEWED_NAME).read_text())
        self.assertEqual(identity["analysis"]["name"],injected)

    def test_existing_target_traversal_and_symlink_fail_without_partial_batch(self):
        reviewed=self._confirm(); reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed),encoding="utf-8")
        (self.root/QMD_NAME).write_text("SYNTHETIC existing",encoding="utf-8")
        with self.assertRaises(IntakeError): generate_scaffold(self.root,reviewed_path)
        self.assertEqual(list(self.root.glob("facs_analysis_*.json")),[])
        outside=Path(self.temp.name).parent/"SYNTHETIC_OUTSIDE_REVIEWED.json"
        outside.write_text(json.dumps(reviewed),encoding="utf-8")
        try:
            with self.assertRaises(IntakeError): generate_scaffold(self.root,outside)
            link=self.root/"SYNTHETIC_link.json"; link.symlink_to(outside)
            with self.assertRaises(IntakeError): generate_scaffold(self.root,link)
        finally:
            outside.unlink(missing_ok=True)

    def test_reviewed_config_tampering_is_rejected(self):
        reviewed=self._confirm(); reviewed["sample_mapping"][0]["condition"]="TAMPERED"
        reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed),encoding="utf-8")
        with self.assertRaises(IntakeError): generate_scaffold(self.root,reviewed_path)

    def test_same_path_input_tamper_and_stale_identity_are_rejected(self):
        reviewed=self._confirm(); reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed),encoding="utf-8")
        (self.root/"SYNTHETIC_sample.fcs").write_text("SYNTHETIC changed",encoding="utf-8")
        with self.assertRaises(IntakeError): generate_scaffold(self.root,reviewed_path)
        template=Path(__file__).parents[1]/"templates/reviewed_analysis.qmd.tmpl"
        template_digest=hashlib.sha256(template.read_bytes()).hexdigest()
        identity={"schema_version":"analysis-identity-1.0","reviewed_config":REVIEWED_NAME,
            "reviewed_content_sha256":reviewed["confirmation"]["content_sha256"],
            "input_provenance_sha256":"0"*64,"template_sha256":template_digest,
            "analysis":reviewed["analyses"][0],"channels":reviewed["channels"],
            "sample_mapping":reviewed["sample_mapping"],"control_relationships":reviewed["control_relationships"],
            "method_choices":{"threshold_specification":None,"normalization_method":None,"background_quantile":None,"dna_alignment_method":None},
            "method_authorization":METHOD_AUTHORIZATION}
        with self.assertRaises(IntakeError): validate_identity(identity,reviewed,REVIEWED_NAME,template_digest)

    def test_batch_write_failure_rolls_back_created_targets(self):
        reviewed=self._confirm(); reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed),encoding="utf-8")
        real_link=os.link; calls=0
        def fail_second(source,target):
            nonlocal calls
            calls += 1
            if calls == 2: raise OSError("SYNTHETIC injected write failure")
            return real_link(source,target)
        with patch("local_facs_assistant.reviewed.os.link",side_effect=fail_second):
            with self.assertRaises(OSError): generate_scaffold(self.root,reviewed_path)
        self.assertFalse((self.root/QMD_NAME).exists())
        self.assertEqual(list(self.root.glob("facs_analysis_*.json")),[])

    def test_rollback_does_not_delete_replaced_target(self):
        reviewed=self._confirm(); reviewed_path=self.root/REVIEWED_NAME
        reviewed_path.write_text(json.dumps(reviewed),encoding="utf-8")
        real_link=os.link; first_target=None; calls=0
        def replace_then_fail(source,target):
            nonlocal calls,first_target
            calls += 1
            if calls == 1:
                first_target=Path(target); return real_link(source,target)
            first_target.unlink(); first_target.write_text("SYNTHETIC replacement",encoding="utf-8")
            raise OSError("SYNTHETIC injected failure after replacement")
        with patch("local_facs_assistant.reviewed.os.link",side_effect=replace_then_fail):
            with self.assertRaises(OSError): generate_scaffold(self.root,reviewed_path)
        self.assertEqual(first_target.read_text(encoding="utf-8"),"SYNTHETIC replacement")


if __name__ == "__main__": unittest.main()
