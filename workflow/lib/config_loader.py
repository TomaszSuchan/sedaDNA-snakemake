"""
Configuration Loader for sedaDNA-snakemake

Simple configuration loading with:
- Structure validation (prevents crashes)
- Required field checking (prevents crashes)
- Helpful error messages (improves UX)
"""

import yaml
from pathlib import Path
from typing import Dict
from copy import deepcopy


# Configuration keys (eliminates magic strings)
class ConfigKeys:
    """Configuration file structure keys"""
    GLOBAL = "global"
    PROJECTS = "projects"
    REFERENCE_DBS = "reference_dbs"
    PARAMETERS = "parameters"
    BARCODES = "barcodes"
    LIBRARIES = "libraries"
    META_ANALYSIS = "meta_analysis"
    FORWARD = "forward"
    REVERSE = "reverse"
    BARCODE_FILE = "barcode_file"


class ConfigurationError(Exception):
    """Raised when configuration is invalid"""
    pass


def deep_merge_dicts(base: Dict, override: Dict) -> Dict:
    """
    Recursively merge two dictionaries.

    Args:
        base: Base dictionary
        override: Dictionary with override values

    Returns:
        Merged dictionary (new object, doesn't modify inputs)
    """
    result = deepcopy(base)

    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge_dicts(result[key], value)
        else:
            result[key] = deepcopy(value)

    return result


def validate_config_structure(config: Dict) -> None:
    """
    Validate basic structure of projects.yaml.

    Args:
        config: Loaded configuration dictionary

    Raises:
        ConfigurationError: If structure is invalid
    """
    # Check required top-level keys
    required_keys = [ConfigKeys.GLOBAL, ConfigKeys.PROJECTS]
    for key in required_keys:
        if key not in config:
            raise ConfigurationError(
                f"Missing required section '{key}' in projects.yaml.\n"
                f"Your configuration must include:\n"
                f"  {ConfigKeys.GLOBAL}:\n"
                f"    {ConfigKeys.REFERENCE_DBS}: ...\n"
                f"    {ConfigKeys.PARAMETERS}: ...\n"
                f"  {ConfigKeys.PROJECTS}:\n"
                f"    ProjectName:\n"
                f"      {ConfigKeys.LIBRARIES}: ..."
            )

    # Validate global section
    global_config = config[ConfigKeys.GLOBAL]
    if not isinstance(global_config, dict):
        raise ConfigurationError(f"'{ConfigKeys.GLOBAL}' must be a dictionary")

    if ConfigKeys.REFERENCE_DBS not in global_config:
        raise ConfigurationError(
            f"Missing '{ConfigKeys.REFERENCE_DBS}' in '{ConfigKeys.GLOBAL}' section"
        )

    if ConfigKeys.PARAMETERS not in global_config:
        raise ConfigurationError(
            f"Missing '{ConfigKeys.PARAMETERS}' in '{ConfigKeys.GLOBAL}' section"
        )

    # Validate projects section
    projects = config[ConfigKeys.PROJECTS]
    if not isinstance(projects, dict):
        raise ConfigurationError(f"'{ConfigKeys.PROJECTS}' must be a dictionary")

    if not projects:
        raise ConfigurationError(
            f"No projects defined. Add at least one project to '{ConfigKeys.PROJECTS}' section."
        )

    # Validate each project
    for proj_name, proj_config in projects.items():
        if not isinstance(proj_config, dict):
            raise ConfigurationError(f"Project '{proj_name}' configuration must be a dictionary")

        if ConfigKeys.LIBRARIES not in proj_config:
            raise ConfigurationError(f"Project '{proj_name}' missing '{ConfigKeys.LIBRARIES}' section")

        libraries = proj_config[ConfigKeys.LIBRARIES]
        if not isinstance(libraries, dict):
            raise ConfigurationError(f"'{ConfigKeys.LIBRARIES}' in project '{proj_name}' must be a dictionary")

        if not libraries:
            raise ConfigurationError(f"Project '{proj_name}' has no libraries defined")

        # Validate each library
        for lib_name, lib_config in libraries.items():
            if not isinstance(lib_config, dict):
                raise ConfigurationError(
                    f"Library '{lib_name}' in project '{proj_name}' must be a dictionary"
                )

            # Check required library keys
            required_lib_keys = [ConfigKeys.FORWARD, ConfigKeys.REVERSE, ConfigKeys.BARCODE_FILE]
            for key in required_lib_keys:
                if key not in lib_config:
                    raise ConfigurationError(
                        f"Library '{lib_name}' in project '{proj_name}' missing required key '{key}'"
                    )


def load_and_validate_config(config_path: str) -> Dict:
    """
    Load and validate projects.yaml configuration.

    Args:
        config_path: Path to projects.yaml

    Returns:
        Validated configuration dictionary

    Raises:
        ConfigurationError: If configuration is invalid
    """
    # Check config file exists
    if not Path(config_path).exists():
        raise ConfigurationError(
            f"Configuration file not found: {config_path}\n"
            f"Please create projects.yaml or use: cp config/projects.example.yaml config/projects.yaml"
        )

    # Load YAML
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Invalid YAML in {config_path}: {e}")

    if config is None:
        raise ConfigurationError(f"Configuration file {config_path} is empty")

    # Validate structure
    validate_config_structure(config)

    return config


def merge_project_parameters(global_params: Dict, project_params: Dict) -> Dict:
    """
    Merge project-specific parameters with global defaults.

    Args:
        global_params: Global parameter defaults
        project_params: Project-specific overrides

    Returns:
        Merged parameters
    """
    if not project_params:
        return deepcopy(global_params)

    return deep_merge_dicts(global_params, project_params)


def extract_projects_info(config: Dict) -> tuple:
    """
    Extract project and library information from config.

    Args:
        config: Validated configuration dictionary

    Returns:
        Tuple of (projects_list, project_libraries_dict, projects_data_dict)
    """
    projects = list(config[ConfigKeys.PROJECTS].keys())

    project_libraries = {
        proj: list(config[ConfigKeys.PROJECTS][proj][ConfigKeys.LIBRARIES].keys())
        for proj in projects
    }

    # Merge parameters for each project
    projects_data = deepcopy(config[ConfigKeys.PROJECTS])
    global_params = config[ConfigKeys.GLOBAL][ConfigKeys.PARAMETERS]

    for proj in projects:
        proj_params = projects_data[proj].get(ConfigKeys.PARAMETERS, {})
        projects_data[proj]["_merged_params"] = merge_project_parameters(
            global_params, proj_params
        )

    return projects, project_libraries, projects_data
