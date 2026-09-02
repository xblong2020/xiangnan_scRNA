import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from scripts.clue_query_module9_8 import (
    map_symbols_to_entrez,
    paired_gmt_text,
    read_gmt_genes,
    safe_job_status,
)


class Module98ClueQueryLogicTest(unittest.TestCase):
    def test_read_and_map_gene_symbols(self):
        with TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input.gmt"
            path.write_text("sig\tdesc\tHNF4A\tSOX4\tMISSING\n", encoding="utf-8")
            genes = read_gmt_genes(path)
        info = pd.DataFrame(
            [
                {"pr_gene_symbol": "HNF4A", "pr_gene_id": 3172, "pr_is_bing": 1},
                {"pr_gene_symbol": "SOX4", "pr_gene_id": 6659, "pr_is_bing": 0},
            ]
        )

        ids, missing, n_bing = map_symbols_to_entrez(genes, info)

        self.assertEqual(ids, ["3172", "6659"])
        self.assertEqual(missing, ["MISSING"])
        self.assertEqual(n_bing, 1)

    def test_paired_gmt_uses_identical_query_name(self):
        text = paired_gmt_text("module9_8", ["3172", "6659"])
        self.assertTrue(text.startswith("module9_8\tmodule9_8_clue_entrez\t"))

    def test_safe_job_status_excludes_secret_and_account_fields(self):
        job = {
            "job_id": "job1",
            "status": "completed",
            "api_key": "secret",
            "user_id": "account",
            "download_url": "private-url",
        }
        safe = safe_job_status(job, Path("result.tar.gz"))

        self.assertNotIn("api_key", safe)
        self.assertNotIn("user_id", safe)
        self.assertNotIn("download_url", safe)
        self.assertFalse(safe["secret_recorded"])


if __name__ == "__main__":
    unittest.main()
