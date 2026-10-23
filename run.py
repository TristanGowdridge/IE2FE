import json
import os
import re

import click
from tqdm import tqdm

from LUSAS_classes import FILENAME_REGEX, LUSASSession

# directory types
read_directory = click.Path(
    exists=True,
    readable=True,
    writable=False,
    file_okay=False,
    dir_okay=True,
    resolve_path=True,
)
write_directory = click.Path(
    exists=True,
    readable=True,
    writable=True,
    file_okay=False,
    dir_okay=True,
    resolve_path=True,
)
# flag types
flag_true = dict(required=False, is_flag=True, flag_value=True, default=None)
flag_false = dict(required=False, is_flag=True, flag_value=False)


@click.command()
@click.option(
    "-ie",
    "--ie-directory",
    prompt="Directory containing IE models to process",
    required=True,
    type=read_directory,
)
@click.option(
    "-fe",
    "--fe-directory",
    prompt="Directory to store the temporary FE models",
    required=True,
    type=write_directory,
)
@click.option(
    "-l",
    "--loadcase-directory",
    prompt="Directory containing the loadcases to run on the models",
    required=True,
    type=read_directory,
)
@click.option(
    "-f",
    "--feature-directory",
    prompt="Directory to save processed features",
    required=True,
    type=write_directory,
)
@click.option(
    "-p",
    "--parameters-file",
    required=False,
    type=click.Path(
        exists=True,
        readable=True,
        writable=False,
        file_okay=True,
        dir_okay=False,
        resolve_path=True,
    ),
)
@click.option("--verbose", required=False, is_flag=True)
@click.option("--limit", required=False, type=int, default=0)
@click.option("--vertical-direction", required=False, type=str)
@click.option("--analysis-category", required=False, type=str)
@click.option("--model-units", required=False, multiple=True, type=str)
@click.option("--use-line-distance", "use_line_distance", **flag_true)
@click.option("--no-use-line-distance", "use_line_distance", **flag_false)
@click.option("--line-mesh-spacing", required=False, multiple=True, type=float)
@click.option("--line-mesh-distance", required=False, multiple=True, type=float)
@click.option("--surface-mesh-spacing", required=False, multiple=True, type=float)
@click.option("--volume-mesh-spacing", required=False, multiple=True, type=float)
@click.option("--n-modal-features", required=False, type=int)
@click.option("--loadcase-format", required=False, type=str)
@click.option("--save-displacement", "save_displacement", **flag_true)
@click.option("--no-save-displacement", "save_displacement", **flag_false)
@click.option("--save-reaction", "save_reaction", **flag_true)
@click.option("--no-save-reaction", "save_reaction", **flag_false)
@click.option("--save-force-moment-beam", "save_force_moment_beam", **flag_true)
@click.option("--no-save-force-moment-beam", "save_force_moment_beam", **flag_false)
@click.option("--save-force-moment-shell", "save_force_moment_shell", **flag_true)
@click.option("--no-save-force-moment-shell", "save_force_moment_shell", **flag_false)
@click.option("--save-loading", "save_loading", **flag_true)
@click.option("--no-save-loading", "save_loading", **flag_false)
@click.option("--save-modal", "save_modal", **flag_true)
@click.option("--no-save-modal", "save_modal", **flag_false)
def run(
    ie_directory: str,
    fe_directory: str,
    loadcase_directory: str,
    feature_directory: str,
    parameters_file: str,
    verbose: bool,
    limit: int,
    **options,
) -> None:

    # print folder settings
    if verbose:
        print(
            "Folder locations:"
            f"\n  IE: {ie_directory}\n  FE: {fe_directory}"
            f"\n  Loadcase: {loadcase_directory}\n  Features: {feature_directory}"
        )

    # find all IE models to process
    filename_re = re.compile(FILENAME_REGEX)
    ie_models = []
    for file in os.listdir(ie_directory):
        if filename_re.fullmatch(file):
            ie_models.append(file)
            if limit > 0 and len(ie_models) >= limit:
                break

    # output number of model diagnostics
    if verbose:
        print(f"Found {len(ie_models)} potential IE models to process")

    # load parameters from JSON file
    parameters = {}
    if parameters_file:
        with open(parameters_file, "r") as fs:
            parameters = json.load(fs)

    # overwrite parameters with any options passed
    for name, value in options.items():
        is_tuple = isinstance(value, tuple)
        if len(name) > 0 and (len(value) > 0 if is_tuple else value != None):
            parts = name.split("_")
            parameters[parts[0] + "".join([p.capitalize() for p in parts[1:]])] = (
                list(value) if is_tuple else value
            )

    # output loaded parameters diagnostics
    if verbose:
        print("FE Parameters:")
        for name, value in parameters.items():
            print(f"  {name}: {value}")

    # process ie models
    LUSASSession.IE_FOLDER = ie_directory
    LUSASSession.FE_FOLDER = fe_directory
    LUSASSession.LOADCASE_FOLDER = loadcase_directory
    LUSASSession.FEATURES_FOLDER = feature_directory
    LUSASSession.VERBOSE = verbose
    LUSASSession.LOGGER = staticmethod(tqdm.write)
    for model in tqdm(ie_models, desc="IE Models"):
        session = LUSASSession.from_ie_model(model)
        session.assign_loadcase()
        session.save_results()
        session.save_model()


if __name__ == "__main__":
    run()
