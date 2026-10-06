"""Reject control-byte and math-markup regressions in public scientific prose."""
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class PublicDocsMarkupTests(unittest.TestCase):
    def test_no_control_bytes_in_prose(self):
        for path in DOCS.rglob("*.md"):
            in_fence = False
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                stripped = line.lstrip()
                if stripped.startswith((chr(96) * 3, "~~~")):
                    in_fence = not in_fence
                    continue
                with self.subTest(path=path.name, line=number):
                    self.assertFalse(
                        any(ord(char) < 32 and char != "\t" for char in line),
                        "unexpected C0 control character",
                    )
                    if not in_fence:
                        self.assertNotIn("\t", line.lstrip(" \t"),
                                         "tab inside prose")

    def test_corrected_equations_have_math_delimiters(self):
        expected = {
            "foundations.md": (
                r"$\operatorname{Tr}(X)I/4$",
                r"$\tau_{KM}$",
                r"$S=\begin{pmatrix}a&b\\b&-a\end{pmatrix}$",
            ),
            "quark-processes.md": (
                r"$g=\operatorname{diag}(1,-1,-1,-1)$",
                r"$\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3$",
                r"$\gamma^{\pm}=(\gamma^0\pm\gamma^3)/\sqrt{2}$",
                r"$\bar{\mathcal G}_B=-\lambda_B g_{1L}^{\bar q/B}$",
                r"$\rho^\Theta=U_T\rho^*U_T^\dagger$",
            ),
        }
        for name, expressions in expected.items():
            source = (DOCS / name).read_text(encoding="utf-8")
            for expression in expressions:
                with self.subTest(path=name, expression=expression):
                    self.assertIn(expression, source)


if __name__ == "__main__":
    unittest.main()
