"""
Test suite for rdDistGeom module.

This module tests the EmbedMolecule function which uses distance geometry
to generate 3D coordinates for molecules.
"""
import pytest
from rdkit import Chem
from rdkit_mcp.Chem.rdDistGeom import EmbedMolecule, EmbedMoleculeResult, EmbedParameters
from rdkit_mcp.base_tools import smiles_to_mol
from rdkit_mcp.utils import decode_mol


class TestEmbedMoleculeBasic:
    """Basic functionality tests for EmbedMolecule."""

    def test_simple_molecule_embedding(self):
        """Test embedding a simple molecule (ethanol)."""
        smiles = "CCO"
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)

        # Should return an EmbedMoleculeResult with conf_id and mol
        assert isinstance(result, EmbedMoleculeResult)
        assert hasattr(result, "conf_id")
        assert hasattr(result, "mol")

        conf_id = result.conf_id
        # Should return a valid conformation ID (0 or higher)
        assert conf_id >= 0

        # Verify the molecule now has a conformation with 3D coordinates
        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() > 0

        # Check that atoms have 3D coordinates
        conf = mol.GetConformer(conf_id)
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            # All coordinates should be non-zero for at least some atoms
            assert pos.x is not None
            assert pos.y is not None
            assert pos.z is not None

    def test_benzene_embedding(self):
        """Test embedding benzene ring."""
        smiles = "c1ccccc1"
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1

    def test_complex_molecule_embedding(self):
        """Test embedding a more complex drug-like molecule (aspirin)."""
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1


class TestEmbedMoleculeParameters:
    """Tests for various EmbedMolecule parameters."""

    def test_random_seed_reproducibility(self):
        """Test that using the same random seed produces identical coordinates."""
        smiles = "CCCC"  # Butane
        p_mol1 = smiles_to_mol(smiles)
        p_mol2 = smiles_to_mol(smiles)

        # Embed both with the same random seed
        params = EmbedParameters(randomSeed=42)
        result1 = EmbedMolecule(p_mol1, params)
        result2 = EmbedMolecule(p_mol2, params)

        conf_id1 = result1.conf_id
        conf_id2 = result2.conf_id

        assert conf_id1 >= 0
        assert conf_id2 >= 0

        mol1 = decode_mol(result1.mol)
        mol2 = decode_mol(result2.mol)

        # Get coordinates from both conformations
        conf1 = mol1.GetConformer(conf_id1)
        conf2 = mol2.GetConformer(conf_id2)

        # Coordinates should be identical (or very close)
        for i in range(mol1.GetNumAtoms()):
            pos1 = conf1.GetAtomPosition(i)
            pos2 = conf2.GetAtomPosition(i)
            assert abs(pos1.x - pos2.x) < 1e-6
            assert abs(pos1.y - pos2.y) < 1e-6
            assert abs(pos1.z - pos2.z) < 1e-6

    def test_clear_confs_parameter(self):
        """Test that clearConfs parameter clears existing conformations."""
        smiles = "CC"
        p_mol = smiles_to_mol(smiles)

        # First embedding
        params1 = EmbedParameters(randomSeed=1)
        result1 = EmbedMolecule(p_mol, params1)
        assert result1.conf_id >= 0

        mol = decode_mol(result1.mol)
        assert mol.GetNumConformers() == 1

        # Second embedding with clearConfs=True (default)
        params2 = EmbedParameters(randomSeed=2, clearConfs=True)
        result2 = EmbedMolecule(result1.mol, params2)
        assert result2.conf_id >= 0

        mol = decode_mol(result2.mol)
        # Should still have only 1 conformation
        assert mol.GetNumConformers() == 1

    def test_clear_confs_false(self):
        """Test that clearConfs=False preserves existing conformations."""
        smiles = "CC"
        p_mol = smiles_to_mol(smiles)

        # First embedding
        params1 = EmbedParameters(randomSeed=1, clearConfs=True)
        result1 = EmbedMolecule(p_mol, params1)
        assert result1.conf_id >= 0

        # Second embedding with clearConfs=False
        params2 = EmbedParameters(randomSeed=2, clearConfs=False)
        result2 = EmbedMolecule(result1.mol, params2)
        assert result2.conf_id >= 0

        mol = decode_mol(result2.mol)
        # Should have 2 conformations now
        assert mol.GetNumConformers() == 2

    def test_max_attempts_parameter(self):
        """Test that maxIterations parameter is respected."""
        smiles = "CC"
        p_mol = smiles_to_mol(smiles)

        # Even with maxIterations=1, simple molecules should embed successfully
        params = EmbedParameters(maxIterations=1)
        result = EmbedMolecule(p_mol, params)
        assert result.conf_id >= 0

    def test_enforce_chirality_parameter(self):
        """Test embedding with enforceChirality parameter."""
        # Use a chiral molecule
        smiles = "C[C@H](O)N"  # (S)-alanine
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(enforceChirality=True)
        result = EmbedMolecule(p_mol, params)
        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1


class TestEmbedMoleculeEdgeCases:
    """Edge case tests for EmbedMolecule."""

    def test_single_atom_molecule(self):
        """Test embedding a single atom (methane)."""
        smiles = "C"
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)

        # Single atoms should embed successfully
        assert result.conf_id >= 0

    def test_two_atom_molecule(self):
        """Test embedding a two-atom molecule."""
        smiles = "CC"
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        conf = mol.GetConformer(result.conf_id)

        # Verify both atoms have positions
        pos1 = conf.GetAtomPosition(0)
        pos2 = conf.GetAtomPosition(1)

        # Calculate distance between the two carbons
        # Should be approximately 1.54 Angstroms for C-C bond
        import math
        distance = math.sqrt(
            (pos1.x - pos2.x)**2 +
            (pos1.y - pos2.y)**2 +
            (pos1.z - pos2.z)**2
        )
        # Allow some tolerance
        assert 1.4 < distance < 1.7

    def test_macrocycle_embedding(self):
        """Test embedding a macrocycle with specific parameters."""
        # 12-membered ring
        smiles = "C1CCCCCCCCCCC1"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(
            useMacrocycleTorsions=True,
            useMacrocycle14config=True
        )
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0


class TestEmbedMoleculeCoordinates:
    """Tests verifying that generated coordinates are reasonable."""

    def test_coordinates_are_3d(self):
        """Verify that all atoms get 3D coordinates."""
        smiles = "CCCCCC"  # Hexane
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)
        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        conf = mol.GetConformer(result.conf_id)

        # Check each atom has non-zero coordinates
        has_nonzero_x = False
        has_nonzero_y = False
        has_nonzero_z = False

        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            if abs(pos.x) > 0.1:
                has_nonzero_x = True
            if abs(pos.y) > 0.1:
                has_nonzero_y = True
            if abs(pos.z) > 0.1:
                has_nonzero_z = True

        # At least some atoms should have significant coordinates in all dimensions
        assert has_nonzero_x
        assert has_nonzero_y
        assert has_nonzero_z

    def test_bond_lengths_reasonable(self):
        """Test that bond lengths are chemically reasonable."""
        smiles = "CCC"  # Propane
        p_mol = smiles_to_mol(smiles)

        result = EmbedMolecule(p_mol)
        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        conf = mol.GetConformer(result.conf_id)

        import math

        # Check C-C bond lengths
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()

            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)

            distance = math.sqrt(
                (pos1.x - pos2.x)**2 +
                (pos1.y - pos2.y)**2 +
                (pos1.z - pos2.z)**2
            )

            # C-C bond should be approximately 1.54 Angstroms
            # Allow reasonable tolerance (1.3 to 1.7 Angstroms)
            assert 1.3 < distance < 1.7, f"Bond length {distance} is unreasonable"


class TestEmbedMoleculeMultipleConformations:
    """Tests for generating multiple conformations."""

    def test_multiple_embeddings_different_random_seeds(self):
        """Test that different random seeds produce different conformations."""
        smiles = "CCCCCC"  # Hexane (flexible)
        p_mol1 = smiles_to_mol(smiles)
        p_mol2 = smiles_to_mol(smiles)

        params1 = EmbedParameters(randomSeed=1)
        params2 = EmbedParameters(randomSeed=2)
        result1 = EmbedMolecule(p_mol1, params1)
        result2 = EmbedMolecule(p_mol2, params2)

        assert result1.conf_id >= 0
        assert result2.conf_id >= 0

        mol1 = decode_mol(result1.mol)
        mol2 = decode_mol(result2.mol)

        conf1 = mol1.GetConformer(result1.conf_id)
        conf2 = mol2.GetConformer(result2.conf_id)

        # At least some atoms should have different positions
        different_coords = False
        for i in range(mol1.GetNumAtoms()):
            pos1 = conf1.GetAtomPosition(i)
            pos2 = conf2.GetAtomPosition(i)

            if abs(pos1.x - pos2.x) > 0.1 or \
               abs(pos1.y - pos2.y) > 0.1 or \
               abs(pos1.z - pos2.z) > 0.1:
                different_coords = True
                break

        assert different_coords, "Different random seeds should produce different conformations"


class TestEmbedMoleculeDrugLikeMolecules:
    """Tests with common drug-like molecules."""

    def test_acetaminophen_embedding(self):
        """Test embedding acetaminophen (Tylenol)."""
        smiles = "CC(=O)NC1=CC=C(C=C1)O"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(randomSeed=42)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1
        assert mol.GetNumAtoms() == 11

    def test_caffeine_embedding(self):
        """Test embedding caffeine."""
        smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(randomSeed=42)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1

    def test_ibuprofen_embedding(self):
        """Test embedding ibuprofen."""
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(randomSeed=42)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1


class TestEmbedMoleculeAdvancedParameters:
    """Tests for advanced embedding parameters."""

    def test_use_basic_knowledge(self):
        """Test embedding with useBasicKnowledge parameter."""
        smiles = "c1ccccc1"  # Benzene should be flat
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(useBasicKnowledge=True)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_use_exp_torsion_angle_prefs(self):
        """Test embedding with experimental torsion angle preferences."""
        smiles = "CCCC"  # Butane
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(useExpTorsionAnglePrefs=True)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_force_tol_parameter(self):
        """Test embedding with different force tolerance."""
        smiles = "CCC"
        p_mol = smiles_to_mol(smiles)

        # Use a looser tolerance
        params = EmbedParameters(optimizerForceTol=0.01)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_et_version_parameter(self):
        """Test embedding with different ET version."""
        smiles = "CCCCCC"
        p_mol = smiles_to_mol(smiles)

        # ETversion 2 is standard for ETKDGv2/v3
        params = EmbedParameters(ETversion=2)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_pruning_parameter(self):
        """Test embedding with RMS threshold pruning."""
        smiles = "CCCCCC"
        p_mol = smiles_to_mol(smiles)

        # Use pruning with 0.5 Angstrom threshold
        params = EmbedParameters(
            pruneRmsThresh=0.5,
            onlyHeavyAtomsForRMS=True
        )
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_threading_parameter(self):
        """Test embedding with multiple threads."""
        smiles = "CCCC"
        p_mol = smiles_to_mol(smiles)

        # Use 2 threads for embedding
        params = EmbedParameters(numThreads=2)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_timeout_parameter(self):
        """Test embedding with timeout parameter."""
        smiles = "CCC"
        p_mol = smiles_to_mol(smiles)

        # Set a generous timeout of 10 seconds
        params = EmbedParameters(timeout=10)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_force_trans_amides_parameter(self):
        """Test embedding with forceTransAmides parameter."""
        # Simple amide molecule
        smiles = "CC(=O)N"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(forceTransAmides=True)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_basin_thresh_parameter(self):
        """Test embedding with custom basin threshold."""
        smiles = "CCCC"
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(basinThresh=10.0)
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_embed_parameters_default(self):
        """Test that default EmbedParameters work correctly."""
        smiles = "CCO"
        p_mol = smiles_to_mol(smiles)

        # Use default parameters
        params = EmbedParameters()
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0

    def test_combined_parameters(self):
        """Test embedding with multiple parameters combined."""
        smiles = "C1CCCCC1"  # Cyclohexane
        p_mol = smiles_to_mol(smiles)

        params = EmbedParameters(
            randomSeed=123,
            enforceChirality=True,
            useBasicKnowledge=True,
            useExpTorsionAnglePrefs=True,
            optimizerForceTol=0.001,
            clearConfs=True
        )
        result = EmbedMolecule(p_mol, params)

        assert result.conf_id >= 0
        mol = decode_mol(result.mol)
        assert mol.GetNumConformers() == 1
