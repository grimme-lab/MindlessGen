*This `CONTRIBUTING` guideline was thankfully taken over from the [BayBE](https://github.com/emdgroup/baybe) code base.*

# Contributing to MindlessGen

**All contributions to MindlessGen are welcome!**

... no matter if bug fixes, new features, or just typo corrections.

To shorten the overall development and review process, this page contains are a
few sections that can make your life easier.

## General Workflow

To implement your contributions in a local development environment,
we recommend the following workflow:

1. Clone a [fork](https://github.com/grimme-lab/MindlessGen/fork) of the repository to
   your local machine.

1. Create and activate a virtual python environment using one of the supported
   python versions.

1. Change into the root folder of the cloned repository and install an editable version
   including all development dependencies:
   ```console
   pip install -e '.[dev]'
   ```

1. Run our tests to verify everything works as expected:
   ```console
   pytest -vv --optional
   ```

1. Install our [pre-commit](https://pre-commit.com/) hooks:
   ```console
   pre-commit install
   ```

1. Create a new branch for your contribution:
   ```console
   git checkout -b <your_branch_name>
   ```

1. **Implement your changes.**

1. Optional but recommended to prevent complaints from our CI pipeline:
   **Test your code.**

   There are several test environments you can run via `tox`, each corresponding to a
   [developer tool](#developer-tools) in a certain Python version.
   You can retrieve all available environments via `tox list`.

   If you want to challenge your machine, you can run all checks in all Python versions
   in parallel via:
   ```console
   tox -p
   ```

   This can be considered the ultimate one-stop check to make sure your code is ready
   for merge.

1. Before setting up a pull request for your contribution, include all of the changes in the `CHANGELOG.md` [file](https://github.com/grimme-lab/MindlessGen/blob/main/CHANGELOG.md).

1. Push the updated branch back to your fork:
   ```console
   git push -u origin <your_branch_name>
   ```

1. Open a pull request via Github's web page.

## Developer Tools

In order to maintain a high code quality, we use a variety of code developer tools.
When following the above described workflow, [pre-commit](https://pre-commit.com/)
will automatically trigger some of these checks during your development process.
In any case, these checks are also conducted in our CI pipeline, which must pass
before your pull request is considered ready for review.
If you have questions or problems, simply ask for advice.

| Tool                                                                                            | Purpose                                   |
|:------------------------------------------------------------------------------------------------|:------------------------------------------|
| [ruff](https://docs.astral.sh/ruff/)                                                            | code linting and formatting               |
| [mypy](https://mypy.readthedocs.io/)                                                            | static type checking                      |
| [pytest](https://docs.pytest.org/)                                                              | testing                                   |
| [tox](https://tox.wiki/)                                                                        | orchestrating all the above               |
| [coverage](https://pypi.org/project/coverage/)                                                  | coverage check and reports                |
