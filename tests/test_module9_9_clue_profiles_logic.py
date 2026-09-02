import unittest

from scripts.clue_landmark_profiles_module9_9 import build_clue_profile_payload


class Module99ClueProfilesLogicTest(unittest.TestCase):
    def test_combined_payload_contains_paired_up_and_down_sets(self):
        payload = build_clue_profile_payload(
            "combined_balanced",
            up_entrez=["3172", "5465"],
            down_entrez=["6659", "3725"],
        )

        self.assertEqual(payload["es_tail"], "both")
        self.assertIn("uptag-cmapfile", payload)
        self.assertIn("dntag-cmapfile", payload)
        self.assertTrue(
            payload["uptag-cmapfile"].startswith(
                "module9_9_combined_balanced\tmodule9_9_clue_entrez\t"
            )
        )
        self.assertTrue(
            payload["dntag-cmapfile"].startswith(
                "module9_9_combined_balanced\tmodule9_9_clue_entrez\t"
            )
        )

    def test_one_sided_payload_omits_empty_gene_set(self):
        malignant = build_clue_profile_payload(
            "malignant_only",
            up_entrez=[],
            down_entrez=["6659", "3725"],
        )
        rescue = build_clue_profile_payload(
            "rescue_only",
            up_entrez=["3172", "5465"],
            down_entrez=[],
        )

        self.assertEqual(malignant["es_tail"], "down")
        self.assertNotIn("uptag-cmapfile", malignant)
        self.assertIn("dntag-cmapfile", malignant)
        self.assertEqual(rescue["es_tail"], "up")
        self.assertIn("uptag-cmapfile", rescue)
        self.assertNotIn("dntag-cmapfile", rescue)


if __name__ == "__main__":
    unittest.main()
