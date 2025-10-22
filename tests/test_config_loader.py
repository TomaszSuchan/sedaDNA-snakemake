"""
Unit tests for config_loader module

Tests configuration loading, validation, and security features.
"""

import pytest
import tempfile
import os
from pathlib import Path
import sys

# Add workflow/lib to path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow/lib"))

from config_loader import (
    validate_project_name,
    validate_library_name,
    validate_file_path,
    validate_file_size,
    deep_merge_dicts,
    validate_config_structure,
    load_and_validate_config,
    ConfigurationError,
    SecurityError,
    ConfigKeys
)


class TestProjectNameValidation:
    """Test project name validation and sanitization"""

    def test_valid_project_names(self):
        """Valid project names should pass"""
        valid_names = ["Project1", "Site_A", "Alps-2023", "Test_123"]
        for name in valid_names:
            assert validate_project_name(name) == name

    def test_invalid_characters(self):
        """Project names with invalid characters should be rejected"""
        invalid_names = [
            "Project;DROP TABLE",  # SQL injection attempt
            "Project../etc",       # Path traversal
            "Project<script>",     # XSS attempt
            "Project name",        # Space
            "Project!",            # Special char
        ]
        for name in invalid_names:
            with pytest.raises(SecurityError, match="Invalid project name"):
                validate_project_name(name)

    def test_starts_with_hyphen(self):
        """Project names starting with hyphen should be rejected"""
        with pytest.raises(SecurityError, match="cannot start with hyphen"):
            validate_project_name("-Project")

    def test_empty_name(self):
        """Empty project name should be rejected"""
        with pytest.raises(SecurityError, match="non-empty string"):
            validate_project_name("")

    def test_too_long(self):
        """Project names over 100 characters should be rejected"""
        long_name = "A" * 101
        with pytest.raises(SecurityError, match="too long"):
            validate_project_name(long_name)

    def test_non_string(self):
        """Non-string project names should be rejected"""
        with pytest.raises(SecurityError):
            validate_project_name(123)


class TestLibraryNameValidation:
    """Test library name validation"""

    def test_valid_library_names(self):
        """Valid library names should pass"""
        valid_names = ["LIB1", "Sample_A", "MET-1"]
        for name in valid_names:
            assert validate_library_name(name) == name

    def test_invalid_characters(self):
        """Library names with invalid characters should be rejected"""
        with pytest.raises(SecurityError, match="Invalid library name"):
            validate_library_name("lib/etc/passwd")


class TestFilePathValidation:
    """Test file path validation and security"""

    def test_valid_absolute_path(self, tmp_path):
        """Valid absolute paths should be accepted"""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        validated = validate_file_path(str(test_file), must_exist=True)
        assert Path(validated).exists()

    def test_file_not_exists_with_check(self):
        """Non-existent file with must_exist=True should raise error"""
        with pytest.raises(FileNotFoundError):
            validate_file_path("/nonexistent/file.txt", must_exist=True)

    def test_directory_not_file(self, tmp_path):
        """Directory paths should be rejected when file expected"""
        with pytest.raises(SecurityError, match="not a file"):
            validate_file_path(str(tmp_path), must_exist=True)

    def test_empty_path(self):
        """Empty path should be rejected"""
        with pytest.raises(SecurityError, match="non-empty string"):
            validate_file_path("")


class TestFileSizeValidation:
    """Test file size limits"""

    def test_large_file_rejected(self, tmp_path):
        """Files exceeding size limit should be rejected"""
        large_file = tmp_path / "large.txt"
        # Create 10MB file
        with open(large_file, 'wb') as f:
            f.write(b'0' * (10 * 1024 * 1024 + 1))

        with pytest.raises(SecurityError, match="exceeds limit"):
            validate_file_size(str(large_file), max_size_mb=10)

    def test_small_file_accepted(self, tmp_path):
        """Files within size limit should be accepted"""
        small_file = tmp_path / "small.txt"
        small_file.write_text("test")

        validate_file_size(str(small_file), max_size_mb=1)  # Should not raise


class TestDeepMergeDicts:
    """Test deep dictionary merging"""

    def test_simple_merge(self):
        """Simple non-nested merge"""
        base = {"a": 1, "b": 2}
        override = {"b": 3, "c": 4}
        result = deep_merge_dicts(base, override)

        assert result == {"a": 1, "b": 3, "c": 4}
        assert base == {"a": 1, "b": 2}  # Original unchanged

    def test_nested_merge(self):
        """Nested dictionary merge"""
        base = {"level1": {"a": 1, "b": 2}}
        override = {"level1": {"b": 3, "c": 4}}
        result = deep_merge_dicts(base, override)

        assert result == {"level1": {"a": 1, "b": 3, "c": 4}}

    def test_deep_nested_merge(self):
        """Multi-level nested merge"""
        base = {
            "level1": {
                "level2": {"a": 1, "b": 2},
                "other": "value"
            }
        }
        override = {
            "level1": {
                "level2": {"b": 3, "c": 4}
            }
        }
        result = deep_merge_dicts(base, override)

        expected = {
            "level1": {
                "level2": {"a": 1, "b": 3, "c": 4},
                "other": "value"
            }
        }
        assert result == expected


class TestConfigStructureValidation:
    """Test configuration structure validation"""

    def test_missing_global_section(self):
        """Config without 'global' section should fail"""
        config = {"projects": {}}
        with pytest.raises(ConfigurationError, match="Missing required section 'global'"):
            validate_config_structure(config)

    def test_missing_projects_section(self):
        """Config without 'projects' section should fail"""
        config = {"global": {}}
        with pytest.raises(ConfigurationError, match="Missing required section 'projects'"):
            validate_config_structure(config)

    def test_empty_projects(self):
        """Config with no projects should fail"""
        config = {
            "global": {
                "reference_dbs": {},
                "parameters": {}
            },
            "projects": {}
        }
        with pytest.raises(ConfigurationError, match="No projects defined"):
            validate_config_structure(config)

    def test_missing_reference_dbs(self):
        """Global section without reference_dbs should fail"""
        config = {
            "global": {"parameters": {}},
            "projects": {"P1": {}}
        }
        with pytest.raises(ConfigurationError, match="reference_dbs"):
            validate_config_structure(config)

    def test_project_without_libraries(self):
        """Project without libraries section should fail"""
        config = {
            "global": {
                "reference_dbs": {"DB1": "path"},
                "parameters": {}
            },
            "projects": {"Project1": {}}
        }
        with pytest.raises(ConfigurationError, match="missing 'libraries'"):
            validate_config_structure(config)

    def test_empty_libraries(self):
        """Project with no libraries should fail"""
        config = {
            "global": {
                "reference_dbs": {"DB1": "path"},
                "parameters": {}
            },
            "projects": {
                "Project1": {
                    "libraries": {}
                }
            }
        }
        with pytest.raises(ConfigurationError, match="no libraries defined"):
            validate_config_structure(config)

    def test_library_missing_required_keys(self):
        """Library without required keys should fail"""
        config = {
            "global": {
                "reference_dbs": {"DB1": "path"},
                "parameters": {}
            },
            "projects": {
                "Project1": {
                    "libraries": {
                        "Lib1": {"forward": "path"}  # Missing reverse and barcode_file
                    }
                }
            }
        }
        with pytest.raises(ConfigurationError, match="missing required key"):
            validate_config_structure(config)

    def test_valid_structure(self):
        """Valid configuration structure should pass"""
        config = {
            "global": {
                "reference_dbs": {"DB1": "path"},
                "parameters": {"param1": "value"}
            },
            "projects": {
                "Project1": {
                    "libraries": {
                        "Lib1": {
                            "forward": "forward.fq.gz",
                            "reverse": "reverse.fq.gz",
                            "barcode_file": "barcodes.txt"
                        }
                    }
                }
            }
        }
        validate_config_structure(config)  # Should not raise


class TestLoadAndValidateConfig:
    """Test full configuration loading"""

    def test_missing_config_file(self):
        """Non-existent config file should fail"""
        with pytest.raises(ConfigurationError, match="not found"):
            load_and_validate_config("/nonexistent/config.yaml")

    def test_empty_config_file(self, tmp_path):
        """Empty config file should fail"""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("")

        with pytest.raises(ConfigurationError, match="empty"):
            load_and_validate_config(str(config_file))

    def test_invalid_yaml(self, tmp_path):
        """Invalid YAML should fail"""
        config_file = tmp_path / "config.yaml"
        config_file.write_text("invalid: yaml: content: [")

        with pytest.raises(ConfigurationError, match="Invalid YAML"):
            load_and_validate_config(str(config_file))

    def test_yaml_bomb_protection(self, tmp_path):
        """Excessively large config files should be rejected"""
        config_file = tmp_path / "config.yaml"
        # Create 15MB file
        with open(config_file, 'w') as f:
            f.write("a: " + "x" * (15 * 1024 * 1024))

        with pytest.raises(SecurityError, match="Maximum size"):
            load_and_validate_config(str(config_file))

    def test_valid_config(self, tmp_path):
        """Valid configuration should load successfully"""
        config_file = tmp_path / "config.yaml"
        config_content = """
global:
  reference_dbs:
    DB1: {db_path}
  parameters:
    max-cpu: 4
  barcodes:
    matching: strict

projects:
  TestProject:
    libraries:
      TestLib:
        forward: {fwd_path}
        reverse: {rev_path}
        barcode_file: {bc_path}
""".format(
            db_path=tmp_path / "db.fasta",
            fwd_path=tmp_path / "fwd.fq.gz",
            rev_path=tmp_path / "rev.fq.gz",
            bc_path=tmp_path / "barcodes.txt"
        )
        config_file.write_text(config_content)

        # Create dummy files
        (tmp_path / "db.fasta").write_text(">test\\nACGT")
        (tmp_path / "fwd.fq.gz").write_bytes(b"test")
        (tmp_path / "rev.fq.gz").write_bytes(b"test")
        (tmp_path / "barcodes.txt").write_text("test")

        config = load_and_validate_config(str(config_file), check_file_existence=False)

        assert ConfigKeys.GLOBAL in config
        assert ConfigKeys.PROJECTS in config
        assert "TestProject" in config[ConfigKeys.PROJECTS]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
