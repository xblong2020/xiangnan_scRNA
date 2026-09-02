import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from scripts.export_codex_project_sessions_to_markdown import export_project_sessions


class ExportCodexProjectSessionsTest(unittest.TestCase):
    def test_export_keeps_visible_dialogue_and_redacts_local_context(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            project_root = root / "project"
            project_root.mkdir()
            codex_home = root / ".codex"
            session_dir = codex_home / "sessions" / "2026" / "07" / "01"
            session_dir.mkdir(parents=True)
            output_dir = project_root / "docs" / "chatgpt_discussions" / "fixture"
            session_id = "00000000-0000-0000-0000-000000000001"

            records = [
                {
                    "timestamp": "2026-07-01T01:00:00Z",
                    "type": "session_meta",
                    "payload": {"id": session_id, "cwd": str(project_root)},
                },
                {
                    "timestamp": "2026-07-01T01:01:00Z",
                    "type": "event_msg",
                    "payload": {"type": "user_message", "message": "请检查 SOX4 结果。"},
                },
                {
                    "timestamp": "2026-07-01T01:01:30Z",
                    "type": "response_item",
                    "payload": {"type": "message", "role": "developer", "content": "internal-only"},
                },
                {
                    "timestamp": "2026-07-01T01:02:00Z",
                    "type": "event_msg",
                    "payload": {
                        "type": "agent_message",
                        "phase": "final",
                        "message": (
                            f"结果位于 {project_root}；输入 C:/private-input.pdf；"
                            "安装路径C:\\Program Files\\R\\R-4.6.0；"
                            "密钥 sk-abcdefghijklmnopqrstuvwxyz123456。"
                        ),
                    },
                },
            ]
            session_file = session_dir / f"rollout-2026-07-01T01-00-00-{session_id}.jsonl"
            session_file.write_text(
                "\n".join(json.dumps(record, ensure_ascii=False) for record in records) + "\n",
                encoding="utf-8",
            )
            (codex_home / "session_index.jsonl").write_text(
                json.dumps(
                    {
                        "id": session_id,
                        "thread_name": f"SOX4 讨论 {project_root}",
                        "updated_at": "2026-07-01T01:02:00Z",
                    },
                    ensure_ascii=False,
                )
                + "\n",
                encoding="utf-8",
            )

            report = export_project_sessions(
                codex_home=codex_home,
                project_root=project_root,
                output_dir=output_dir,
            )

            self.assertEqual(report["session_count"], 1)
            self.assertEqual(report["sessions"][0]["user_message_count"], 1)
            self.assertEqual(report["sessions"][0]["assistant_message_count"], 1)
            transcript = (output_dir / report["sessions"][0]["output_file"]).read_text(encoding="utf-8")
            self.assertIn("# SOX4 讨论 $PROJECT_ROOT", transcript)
            self.assertIn("请检查 SOX4 结果。", transcript)
            self.assertIn("$PROJECT_ROOT", transcript)
            self.assertIn("<LOCAL_PATH>", transcript)
            self.assertNotIn(str(project_root), transcript)
            self.assertNotIn("C:/private-input.pdf", transcript)
            self.assertNotIn(r"C:\Program Files\R\R-4.6.0", transcript)
            self.assertNotIn("sk-abcdefghijklmnopqrstuvwxyz123456", transcript)
            self.assertNotIn("internal-only", transcript)
            self.assertTrue((output_dir / "index.md").is_file())
            self.assertTrue((output_dir / "export_report.json").is_file())


if __name__ == "__main__":
    unittest.main()
