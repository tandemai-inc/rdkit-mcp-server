"""
Test suite for base_tools module.

This module tests the core conversion functions in base_tools, including
SMILES, SMARTS, SDF, and PDB conversions to/from RDKit Mol objects.
"""
from rdkit_mcp.base_tools import smarts_to_mol, smiles_to_mol, mol_to_smiles
from rdkit_mcp.utils import decode_mol
from rdkit import Chem
import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestSmartsToMol:
    """Test suite for smarts_to_mol function."""

    def test_simple_aromatic_ring(self):
        """Test converting a simple aromatic benzene ring SMARTS."""
        smarts = "c1ccccc1"

        encoded_mol = smarts_to_mol(smarts)

        # Verify we get a base64 encoded string back
        assert isinstance(encoded_mol, str)
        assert len(encoded_mol) > 0

        # Decode and verify it's a valid mol object
        mol = decode_mol(encoded_mol)
        assert mol is not None
        assert mol.GetNumAtoms() == 6

    def test_aliphatic_benzene_pattern(self):
        """Test converting an aliphatic benzene SMARTS pattern."""
        smarts = "C1=CC=CC=C1"

        encoded_mol = smarts_to_mol(smarts)

        mol = decode_mol(encoded_mol)
        assert mol is not None
        assert mol.GetNumAtoms() == 6

    def test_carbonyl_group(self):
        """Test converting a carbonyl group SMARTS pattern."""
        smarts = "[C]=[O]"

        encoded_mol = smarts_to_mol(smarts)

        mol = decode_mol(encoded_mol)
        assert mol is not None
        assert mol.GetNumAtoms() == 2

    def test_amide_substructure(self):
        """Test converting an amide SMARTS pattern."""
        smarts = "[NX3][CX3](=[OX1])"

        encoded_mol = smarts_to_mol(smarts)

        mol = decode_mol(encoded_mol)
        assert mol is not None
        # Amide has N-C=O (3 atoms)
        assert mol.GetNumAtoms() == 3

    def test_hydroxyl_group(self):
        """Test converting a hydroxyl group SMARTS pattern."""
        smarts = "[OX2H]"

        encoded_mol = smarts_to_mol(smarts)

        mol = decode_mol(encoded_mol)
        assert mol is not None
        assert mol.GetNumAtoms() == 1

    def test_invalid_smarts_raises_error(self):
        """Test that invalid SMARTS raises ToolError."""
        from mcp.server.fastmcp.exceptions import ToolError

        invalid_smarts = "INVALID[[[SMARTS"

        with pytest.raises(ToolError) as exc_info:
            smarts_to_mol(invalid_smarts)

        assert "Invalid or unparsable SMARTS" in str(exc_info.value)

    def test_complex_smarts_pattern(self):
        """Test converting a complex SMARTS pattern with multiple features."""
        # Aromatic carbon with at least one hydrogen
        smarts = "[cH1]"

        encoded_mol = smarts_to_mol(smarts)

        mol = decode_mol(encoded_mol)
        assert mol is not None
        assert mol.GetNumAtoms() == 1


class TestSmartsToMolIntegration:
    """Integration tests for smarts_to_mol with substructure matching."""

    def test_smarts_matches_molecule(self):
        """Test that a SMARTS pattern can be used to match a molecule."""
        # Convert SMARTS to mol
        smarts = "c1ccccc1"  # Aromatic benzene
        smarts_mol_encoded = smarts_to_mol(smarts)
        smarts_mol = decode_mol(smarts_mol_encoded)

        # Create a benzene molecule
        smiles = "c1ccccc1"
        smiles_mol_encoded = smiles_to_mol(smiles)
        smiles_mol = decode_mol(smiles_mol_encoded)

        # Test substructure match
        match = smiles_mol.HasSubstructMatch(smarts_mol)
        assert match is True

    def test_smarts_no_match_different_molecule(self):
        """Test that a SMARTS pattern doesn't match an unrelated molecule."""
        # Aromatic ring pattern
        smarts = "c1ccccc1"
        smarts_mol_encoded = smarts_to_mol(smarts)
        smarts_mol = decode_mol(smarts_mol_encoded)

        # Ethanol (no aromatic ring)
        smiles = "CCO"
        smiles_mol_encoded = smiles_to_mol(smiles)
        smiles_mol = decode_mol(smiles_mol_encoded)

        # Should not match
        match = smiles_mol.HasSubstructMatch(smarts_mol)
        assert match is False

    def test_smarts_matches_substructure(self):
        """Test SMARTS pattern matches as substructure in larger molecule."""
        # Carbonyl pattern
        smarts = "C=O"
        smarts_mol_encoded = smarts_to_mol(smarts)
        smarts_mol = decode_mol(smarts_mol_encoded)

        # Acetone (has carbonyl)
        smiles = "CC(=O)C"
        smiles_mol_encoded = smiles_to_mol(smiles)
        smiles_mol = decode_mol(smiles_mol_encoded)

        # Should match
        match = smiles_mol.HasSubstructMatch(smarts_mol)
        assert match is True

        # Get the matching atoms
        match_atoms = smiles_mol.GetSubstructMatch(smarts_mol)
        assert len(match_atoms) == 2  # C and O


class TestSmartsRoundTrip:
    """Test round-trip conversions involving SMARTS."""

    def test_smarts_to_mol_preserves_pattern(self):
        """Test that SMARTS pattern is preserved through conversion."""
        smarts_patterns = [
            "c1ccccc1",  # Aromatic benzene
            "[NX3][CX3](=[OX1])",  # Amide
            "[OH]",  # Hydroxyl
            "C=O",  # Carbonyl
        ]

        for smarts in smarts_patterns:
            encoded_mol = smarts_to_mol(smarts)
            mol = decode_mol(encoded_mol)

            # Verify the mol object can be created from the original SMARTS
            expected_mol = Chem.MolFromSmarts(smarts)

            # Both should have the same number of atoms
            assert mol.GetNumAtoms() == expected_mol.GetNumAtoms()
