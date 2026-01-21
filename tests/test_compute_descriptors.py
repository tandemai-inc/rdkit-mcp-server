"""
Test dataset for compute_descriptors function.

This module tests the compute_descriptors function which computes
RDKit molecular property descriptors for a collection of ligands using
the efficient rdMolDescriptors.Properties API.
"""
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from rdkit_mcp.Chem.rdMolDescriptors import compute_descriptors


# Test molecules with known descriptor values
# Acetaminophen (Tylenol): CC(=O)NC1=CC=C(C=C1)O
# Aspirin: CC(=O)OC1=CC=CC=C1C(=O)O
# Ethanol: CCO
# Benzene: c1ccccc1
# Caffeine: Cn1cnc2c1c(=O)n(c(=O)n2C)C


class TestComputeDescriptorsRdkitBasic:
    """Basic functionality tests for compute_descriptors."""

    def test_single_molecule_single_descriptor(self):
        """Test computing a single descriptor for a single molecule."""
        smiles_list = ["CCO"]  # Ethanol
        descriptor_names = ["exactmw"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        assert len(result[0]) == 1
        # Ethanol exact MW is approximately 46.042
        assert abs(result[0][0] - 46.042) < 0.01

    def test_single_molecule_multiple_descriptors(self):
        """Test computing multiple descriptors for a single molecule."""
        smiles_list = ["c1ccccc1"]  # Benzene
        descriptor_names = ["exactmw", "NumHeavyAtoms", "NumRotatableBonds"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        assert len(result[0]) == 3
        # Benzene: MW ~78.047, 6 heavy atoms, 0 rotatable bonds
        assert abs(result[0][0] - 78.047) < 0.01
        assert result[0][1] == 6
        assert result[0][2] == 0

    def test_multiple_molecules_single_descriptor(self):
        """Test computing a single descriptor for multiple molecules."""
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]  # Ethanol, Benzene, Acetic acid
        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 3
        # Ethanol: 3 heavy atoms (2 C + 1 O)
        assert result[0][0] == 3
        # Benzene: 6 heavy atoms
        assert result[1][0] == 6
        # Acetic acid: 4 heavy atoms (2 C + 2 O)
        assert result[2][0] == 4

    def test_multiple_molecules_multiple_descriptors(self):
        """Test computing multiple descriptors for multiple molecules."""
        smiles_list = ["CCO", "c1ccccc1"]  # Ethanol, Benzene
        descriptor_names = ["exactmw", "NumHeavyAtoms", "NumAromaticRings"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 2
        assert len(result[0]) == 3
        assert len(result[1]) == 3

        # Ethanol: MW ~46, 3 heavy atoms, 0 aromatic rings
        assert abs(result[0][0] - 46.042) < 0.01
        assert result[0][1] == 3
        assert result[0][2] == 0

        # Benzene: MW ~78, 6 heavy atoms, 1 aromatic ring
        assert abs(result[1][0] - 78.047) < 0.01
        assert result[1][1] == 6
        assert result[1][2] == 1


class TestComputeDescriptorsRdkitEdgeCases:
    """Edge case tests for compute_descriptors."""

    def test_empty_smiles_list(self):
        """Test with an empty SMILES list."""
        smiles_list = []
        descriptor_names = ["exactmw"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert result == []

    def test_invalid_smiles_skipped(self):
        """Test that invalid SMILES are silently skipped."""
        smiles_list = ["CCO", "INVALID_SMILES", "c1ccccc1"]
        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_list, descriptor_names)

        # Only 2 valid molecules should be in the result
        assert len(result) == 2
        assert result[0][0] == 3  # Ethanol
        assert result[1][0] == 6  # Benzene

    def test_all_invalid_smiles(self):
        """Test with all invalid SMILES strings."""
        smiles_list = ["NOT_A_SMILES", "ALSO_INVALID", "12345"]
        descriptor_names = ["exactmw"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert result == []

    def test_generator_input(self):
        """Test that generator input works correctly."""
        def smiles_generator():
            yield "CCO"
            yield "c1ccccc1"

        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_generator(), descriptor_names)

        assert len(result) == 2
        assert result[0][0] == 3
        assert result[1][0] == 6

    def test_tuple_input(self):
        """Test that tuple input works correctly."""
        smiles_tuple = ("CCO", "c1ccccc1")
        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_tuple, descriptor_names)

        assert len(result) == 2

    def test_invalid_descriptor_raises_error(self):
        """Test that invalid descriptor names raise ValueError."""
        smiles_list = ["CCO"]
        descriptor_names = ["NotARealDescriptor"]

        with pytest.raises(Exception):  # RDKit may raise different exception types
            compute_descriptors(smiles_list, descriptor_names)


class TestComputeDescriptorsRdkitDruglikeMolecules:
    """Tests with drug-like molecules to verify descriptor accuracy."""

    def test_acetaminophen_descriptors(self):
        """Test descriptors for acetaminophen (Tylenol)."""
        # CC(=O)NC1=CC=C(C=C1)O - Acetaminophen
        smiles_list = ["CC(=O)NC1=CC=C(C=C1)O"]
        descriptor_names = [
            "exactmw",
            "NumHeavyAtoms",
            "NumHBD",  # H-bond donors
            "NumHBA",  # H-bond acceptors
            "NumRotatableBonds",
            "NumAromaticRings",
        ]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        row = result[0]

        # Acetaminophen: C8H9NO2, MW ~151.063
        assert abs(row[0] - 151.063) < 0.01
        # 11 heavy atoms (8 C + 1 N + 2 O)
        assert row[1] == 11
        # 2 H-bond donors (NH and OH)
        assert row[2] == 2
        # 2 H-bond acceptors (N and O in carbonyl, plus hydroxyl O)
        assert row[3] == 2
        # 1 rotatable bond
        assert row[4] == 1
        # 1 aromatic ring
        assert row[5] == 1

    def test_aspirin_descriptors(self):
        """Test descriptors for aspirin."""
        # CC(=O)OC1=CC=CC=C1C(=O)O - Aspirin
        smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O"]
        descriptor_names = [
            "exactmw",
            "NumHeavyAtoms",
            "NumHBD",
            "NumHBA",
            "NumRotatableBonds",
        ]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        row = result[0]

        # Aspirin: C9H8O4, MW ~180.042
        assert abs(row[0] - 180.042) < 0.01
        # 13 heavy atoms
        assert row[1] == 13
        # 1 H-bond donor (carboxylic acid OH)
        assert row[2] == 1
        # 3 H-bond acceptors
        assert row[3] == 3

    def test_caffeine_descriptors(self):
        """Test descriptors for caffeine."""
        # Cn1cnc2c1c(=O)n(c(=O)n2C)C - Caffeine
        smiles_list = ["Cn1cnc2c1c(=O)n(c(=O)n2C)C"]
        descriptor_names = [
            "exactmw",
            "NumHeavyAtoms",
            "NumHBD",
            "NumHBA",
            "NumRings",
        ]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        row = result[0]

        # Caffeine: C8H10N4O2, MW ~194.080
        assert abs(row[0] - 194.080) < 0.01
        # 14 heavy atoms
        assert row[1] == 14
        # 0 H-bond donors (all N are methylated)
        assert row[2] == 0
        # 2 aromatic + 1 fused = 2 rings total in RDKit counting
        assert row[4] == 2


class TestComputeDescriptorsRdkitDescriptorCoverage:
    """Tests covering various descriptor types available through Properties."""

    def test_molecular_weight_descriptors(self):
        """Test molecular weight related descriptors."""
        smiles_list = ["CCO"]  # Ethanol
        descriptor_names = ["exactmw", "amw"]  # Exact MW and average MW

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        # Exact MW should be close to 46.042
        assert abs(result[0][0] - 46.042) < 0.01
        # Average MW should be close to 46.07 (including average isotope masses)
        assert abs(result[0][1] - 46.07) < 0.1

    def test_ring_descriptors(self):
        """Test ring-related descriptors."""
        smiles_list = ["c1ccc2ccccc2c1"]  # Naphthalene (2 fused aromatic rings)
        descriptor_names = [
            "NumRings",
            "NumAromaticRings",
            "NumAliphaticRings",
        ]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        # Naphthalene: 2 total rings, 2 aromatic, 0 aliphatic
        assert result[0][0] == 2
        assert result[0][1] == 2
        assert result[0][2] == 0

    def test_lipophilicity_descriptors(self):
        """Test lipophilicity-related descriptors."""
        smiles_list = ["CCCCCCCC"]  # Octane (very lipophilic)
        descriptor_names = ["CrippenClogP"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        # Octane should have high LogP (> 3)
        assert result[0][0] > 3.0

    def test_polar_surface_area(self):
        """Test TPSA (topological polar surface area) descriptor."""
        smiles_list = ["c1ccccc1", "CCO"]  # Benzene (no polar), Ethanol (polar)
        descriptor_names = ["tpsa"]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 2
        # Benzene: TPSA should be 0 (no polar atoms)
        assert result[0][0] == 0.0
        # Ethanol: TPSA should be ~20 (from hydroxyl group)
        assert result[1][0] > 15.0


class TestComputeDescriptorsRdkitOrderPreservation:
    """Tests verifying that order is preserved correctly."""

    def test_output_order_matches_input_smiles_order(self):
        """Test that output rows match input SMILES order for valid molecules."""
        # Use molecules with clearly different MW to verify order
        smiles_list = ["C", "CC", "CCC", "CCCC"]  # Methane, Ethane, Propane, Butane
        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_list, descriptor_names)

        # Verify order is preserved
        assert result[0][0] == 1  # Methane: 1 C
        assert result[1][0] == 2  # Ethane: 2 C
        assert result[2][0] == 3  # Propane: 3 C
        assert result[3][0] == 4  # Butane: 4 C

    def test_descriptor_order_matches_input_descriptor_order(self):
        """Test that descriptor values match the order specified."""
        smiles_list = ["c1ccccc1"]  # Benzene
        # Request descriptors in specific order
        descriptor_names_a = ["NumHeavyAtoms", "NumAromaticRings", "exactmw"]
        descriptor_names_b = ["exactmw", "NumAromaticRings", "NumHeavyAtoms"]

        result_a = compute_descriptors(smiles_list, descriptor_names_a)
        result_b = compute_descriptors(smiles_list, descriptor_names_b)

        # First descriptor in result_a (NumHeavyAtoms=6) should match
        # third descriptor in result_b
        assert result_a[0][0] == result_b[0][2]

        # Third descriptor in result_a (exactmw) should match
        # first descriptor in result_b
        assert result_a[0][2] == result_b[0][0]

    def test_invalid_smiles_do_not_affect_order(self):
        """Test that invalid SMILES don't shift the order of valid results."""
        smiles_list = ["C", "INVALID", "CC", "ALSO_INVALID", "CCC"]
        descriptor_names = ["NumHeavyAtoms"]

        result = compute_descriptors(smiles_list, descriptor_names)

        # Should have 3 results for valid SMILES, in original order
        assert len(result) == 3
        assert result[0][0] == 1  # Methane
        assert result[1][0] == 2  # Ethane
        assert result[2][0] == 3  # Propane


class TestComputeDescriptorsRdkitLargeDataset:
    """Performance-oriented tests with larger datasets."""

    def test_batch_of_100_molecules(self):
        """Test processing a batch of 100 molecules."""
        # Generate 100 alkane SMILES
        smiles_list = ["C" * i for i in range(1, 101)]
        descriptor_names = ["exactmw", "NumHeavyAtoms"]

        result = compute_descriptors(smiles_list, descriptor_names)

        # All should be valid
        assert len(result) == 100

        # Verify first and last entries
        # Methane (C): 1 heavy atom
        assert result[0][1] == 1
        # C100: 100 heavy atoms
        assert result[99][1] == 100

    def test_multiple_descriptors_efficiency(self):
        """Test that many descriptors can be computed efficiently."""
        smiles_list = ["CC(=O)NC1=CC=C(C=C1)O"]  # Acetaminophen
        descriptor_names = [
            "exactmw",
            "amw",
            "NumHeavyAtoms",
            "NumHBD",
            "NumHBA",
            "NumRotatableBonds",
            "NumRings",
            "NumAromaticRings",
            "NumAliphaticRings",
            "tpsa",
            "CrippenClogP",
            "CrippenMR",
        ]

        result = compute_descriptors(smiles_list, descriptor_names)

        assert len(result) == 1
        assert len(result[0]) == len(descriptor_names)

        # All values should be valid numbers (not NaN)
        for value in result[0]:
            assert value == value  # NaN check: NaN != NaN
