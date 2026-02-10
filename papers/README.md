# EdgeFEM Technical Papers

This directory contains technical documentation and white papers for the EdgeFEM project.

## Contents

### white_paper/

Technical white paper documenting EdgeFEM's finite element methodology with validated case studies.

**Building the PDF:**

```bash
cd white_paper
make
```

**Requirements:**
- LaTeX distribution (TeX Live, MacTeX, or MiKTeX)
- pdflatex
- bibtex

On macOS:
```bash
brew install --cask mactex-no-gui
# or for minimal installation:
brew install basictex
sudo tlmgr update --self
sudo tlmgr install collection-latexrecommended algorithms siunitx booktabs
```

On Ubuntu/Debian:
```bash
sudo apt-get install texlive-latex-recommended texlive-latex-extra texlive-science
```

**Paper Structure:**
1. Introduction - Edge element advantages, ecosystem positioning
2. Mathematical Formulation - Curl-curl weak form, Nedelec elements
3. Boundary Conditions - PEC, ABC, PML, periodic/Floquet
4. Port Formulation - Eigenvector-based weights, S-parameter extraction
5. Case Study - WR-90 waveguide validation
6. Verification Summary - Test suite coverage
7. Conclusions - Validated capabilities, future work

## Future: IEEE Journal Paper

After the white paper is complete, an IEEE-quality journal paper is planned targeting:
- IEEE Trans. Antennas and Propagation
- IEEE Trans. Microwave Theory and Techniques
- IEEE Antennas and Wireless Propagation Letters

### Additional Work Required for Journal Submission

1. **Extended validation suite**
   - Horn antenna radiation patterns
   - Patch antenna array mutual coupling
   - Frequency-selective surface transmission

2. **Performance benchmarks**
   - Timing vs. mesh size scaling
   - Memory consumption analysis
   - Comparison with commercial tools (if available)

3. **Novel contributions to emphasize**
   - Eigenvector-based port weights (vs. analytical)
   - Optimal ABC scaling factor derivation
   - Python/C++ hybrid architecture for automation

### Timeline (Post White Paper)

- Month 1-2: Extended validation suite
- Month 3: Draft IEEE paper with new results
- Month 4: Co-author review and revision
- Month 5: Submit to target journal

## Citation

If you use EdgeFEM in your research, please cite:

```bibtex
@misc{edgefem2026,
  author = {{VectorEM Project}},
  title = {{EdgeFEM}: A {3D} Finite Element Electromagnetics Simulator},
  year = {2026},
  url = {https://github.com/your-repo/VectorEM}
}
```
