# Contributing to VectorEM

Thank you for your interest in contributing to VectorEM! This document provides guidelines for contributing to the project.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/VectorEM.git
   cd VectorEM
   ```
3. **Build the project** following the [README](README.md) instructions
4. **Run tests** to ensure everything works:
   ```bash
   ctest --test-dir build -j
   ```

## Development Workflow

### Branch Naming

- `feature/description` - New features
- `fix/description` - Bug fixes
- `docs/description` - Documentation updates
- `refactor/description` - Code refactoring

### Making Changes

1. Create a feature branch from `main`:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes, following the code style guidelines below

3. Run formatting and linting:
   ```bash
   bash tools/lint.sh build
   ```

4. Run tests:
   ```bash
   ctest --test-dir build -j
   ```

5. Commit with clear messages:
   ```bash
   git commit -m "Add eigenmode solver for circular waveguides"
   ```

## Code Style

### C++

- **Standard**: C++20
- **Formatting**: clang-format (configuration in `.clang-format`)
- **Static Analysis**: clang-tidy (configuration in `.clang-tidy`)
- **Namespace**: All code in `namespace vectorem`

Run formatting before committing:
```bash
clang-format -i src/*.cpp include/vectorem/*.hpp
```

### Python

- **Version**: Python 3.9+
- **Style**: PEP 8
- **Type Hints**: Encouraged for public APIs

### General Guidelines

- Keep functions focused and under 50 lines when possible
- Use descriptive variable names
- Add comments for non-obvious algorithms
- Prefer `const` correctness in C++

## Testing

### Writing Tests

- Tests use simple `assert()` macros (no external framework)
- Place tests in `tests/` directory
- Name test files `test_*.cpp`
- Add appropriate CTest labels (`smoke`, `pml`, `abc`, etc.)

### Running Tests

```bash
# All tests
ctest --test-dir build -j

# Quick smoke tests
ctest --test-dir build -L smoke

# Specific test
./build/tests/test_maxwell
```

## Pull Request Process

### Before Submitting

1. **Rebase on latest main**:
   ```bash
   git fetch origin
   git rebase origin/main
   ```

2. **Ensure CI will pass**:
   - All tests pass
   - Code is formatted
   - clang-tidy has no warnings

3. **Update documentation** if needed:
   - Update README.md for user-facing changes
   - Update API docs for new functions
   - Add/update tutorials for new features

### PR Requirements

- **Green CI**: All GitHub Actions checks must pass
- **Description**: Clear explanation of what and why
- **Tests**: New features need tests; bug fixes need regression tests
- **Docs**: User-facing changes need documentation updates

### Review Process

1. Submit PR against `main` branch
2. Address reviewer feedback
3. Maintainer merges after approval

## Architecture Decisions

Large features or architectural changes require an **Architecture Decision Record (ADR)**:

1. Create `docs/adr/NNNN-title.md`
2. Describe context, decision, and consequences
3. Get approval before implementation

## Reporting Issues

### Bug Reports

Include:
- VectorEM version (git commit hash)
- Operating system and compiler version
- Minimal reproducing example
- Expected vs actual behavior

### Feature Requests

Include:
- Use case description
- Proposed API (if applicable)
- Alternative approaches considered

## Questions?

- Open a [GitHub Discussion](https://github.com/jman4162/VectorEM/discussions) for questions
- Check existing [Issues](https://github.com/jman4162/VectorEM/issues) for known problems

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
