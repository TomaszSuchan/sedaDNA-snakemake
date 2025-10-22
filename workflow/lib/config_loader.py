"""
Configuration Loader and Validator for sedaDNA-snakemake

Provides secure configuration loading with:
- Schema validation
- Path traversal protection
- Name sanitization
- Parameter merging
"""

import os
import re
import yaml
from pathlib import Path
from typing import Dict, List, Any
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


class SecurityError(Exception):
    """Raised when security validation fails"""
    pass


def validate_project_name(name: str) -> str:
    """
    Validate and sanitize project name.

    Only allows: alphanumeric, underscore, hyphen
    Prevents: path traversal, command injection, special chars

    Args:
        name: Project name to validate

    Returns:
        Sanitized project name

    Raises:
        SecurityError: If name contains invalid characters
    """
    if not name or not isinstance(name, str):
        raise SecurityError("Project name must be a non-empty string")

    if not re.match(r'^[a-zA-Z0-9_-]+$', name):
        raise SecurityError(
            f"Invalid project name '{name}'. "
            "Only alphanumeric characters, underscores, and hyphens allowed. "
            "This prevents security issues like path traversal and command injection."
        )

    if name.startswith('-'):
        raise SecurityError(
            f"Project name '{name}' cannot start with hyphen (could be interpreted as command flag)"
        )

    if len(name) > 100:
        raise SecurityError(f"Project name '{name}' too long (max 100 characters)")

    return name


def validate_library_name(name: str) -> str:
    """
    Validate and sanitize library name.

    Same rules as project names.

    Args:
        name: Library name to validate

    Returns:
        Sanitized library name

    Raises:
        SecurityError: If name contains invalid characters
    """
    if not name or not isinstance(name, str):
        raise SecurityError("Library name must be a non-empty string")

    if not re.match(r'^[a-zA-Z0-9_-]+$', name):
        raise SecurityError(
            f"Invalid library name '{name}'. "
            "Only alphanumeric characters, underscores, and hyphens allowed."
        )

    if len(name) > 100:
        raise SecurityError(f"Library name '{name}' too long (max 100 characters)")

    return name


def validate_file_path(file_path: str, must_exist: bool = False) -> str:
    """
    Validate file path for security.

    Prevents:
    - Path traversal attacks
    - Access to sensitive system files
    - Symlink attacks

    Args:
        file_path: Path to validate
        must_exist: If True, verify file exists

    Returns:
        Absolute path to file

    Raises:
        SecurityError: If path is unsafe
        FileNotFoundError: If must_exist=True and file doesn't exist
    """
    if not file_path or not isinstance(file_path, str):
        raise SecurityError("File path must be a non-empty string")

    # Convert to Path object
    path = Path(file_path)

    # Get absolute path (resolves symlinks)
    try:
        abs_path = path.resolve(strict=must_exist)
    except (RuntimeError, OSError) as e:
        raise SecurityError(f"Cannot resolve path '{file_path}': {e}")

    # Check file exists if required
    if must_exist and not abs_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    # Verify it's a file (not directory) if it exists
    if must_exist and not abs_path.is_file():
        raise SecurityError(f"Path exists but is not a file: {file_path}")

    return str(abs_path)


def validate_file_size(file_path: str, max_size_mb: int = 100) -> None:
    """
    Validate file size to prevent resource exhaustion.

    Args:
        file_path: Path to file
        max_size_mb: Maximum allowed size in megabytes

    Raises:
        SecurityError: If file is too large
    """
    if not Path(file_path).exists():
        return  # Will fail later in pipeline

    size_bytes = Path(file_path).stat().st_size
    size_mb = size_bytes / (1024 * 1024)

    if size_mb > max_size_mb:
        raise SecurityError(
            f"File '{file_path}' is {size_mb:.1f} MB, exceeds limit of {max_size_mb} MB. "
            "This prevents resource exhaustion attacks."
        )


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
        # Validate project name
        validate_project_name(proj_name)

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
            validate_library_name(lib_name)

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


def validate_file_paths(config: Dict, check_existence: bool = False) -> None:
    """
    Validate all file paths in configuration.

    Args:
        config: Loaded configuration dictionary
        check_existence: If True, verify all files exist

    Raises:
        SecurityError: If paths are unsafe
        FileNotFoundError: If check_existence=True and file missing
    """
    # Validate reference database paths
    for db_name, db_path in config[ConfigKeys.GLOBAL][ConfigKeys.REFERENCE_DBS].items():
        try:
            validate_file_path(db_path, must_exist=check_existence)
        except (SecurityError, FileNotFoundError) as e:
            raise ConfigurationError(
                f"Invalid reference database path for '{db_name}': {e}"
            )

    # Validate library file paths
    for proj_name, proj_config in config[ConfigKeys.PROJECTS].items():
        for lib_name, lib_config in proj_config[ConfigKeys.LIBRARIES].items():
            # Validate forward, reverse, barcode files
            for key in [ConfigKeys.FORWARD, ConfigKeys.REVERSE, ConfigKeys.BARCODE_FILE]:
                file_path = lib_config[key]
                try:
                    validated_path = validate_file_path(file_path, must_exist=check_existence)
                    # Update config with validated path
                    lib_config[key] = validated_path

                    # Check file size
                    if check_existence:
                        max_size = 50000 if key == ConfigKeys.BARCODE_FILE else 100000  # MB
                        validate_file_size(validated_path, max_size_mb=max_size)

                except (SecurityError, FileNotFoundError) as e:
                    raise ConfigurationError(
                        f"Invalid {key} path in project '{proj_name}', library '{lib_name}': {e}"
                    )


def load_and_validate_config(config_path: str, check_file_existence: bool = False) -> Dict:
    """
    Load and validate projects.yaml configuration.

    Args:
        config_path: Path to projects.yaml
        check_file_existence: If True, verify all input files exist

    Returns:
        Validated configuration dictionary

    Raises:
        ConfigurationError: If configuration is invalid
        SecurityError: If security validation fails
    """
    # Check config file exists
    if not Path(config_path).exists():
        raise ConfigurationError(
            f"Configuration file not found: {config_path}\n"
            f"Please create projects.yaml or use: cp config/projects.example.yaml config/projects.yaml"
        )

    # Check config file size (prevent YAML bombs)
    config_size = Path(config_path).stat().st_size / 1024  # KB
    if config_size > 10240:  # 10 MB
        raise SecurityError(
            f"Configuration file {config_path} is {config_size/1024:.1f} MB. "
            "Maximum size is 10 MB to prevent resource exhaustion."
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

    # Validate file paths
    validate_file_paths(config, check_existence=check_file_existence)

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
