from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

from nomad.config import config
from nomad.metainfo import (  # pylint: disable=unused-import
    Category,
    MCategory,
    MSection,
    Package,
    Quantity,
    Reference,
    SchemaPackage,
    Section,
    SectionProxy,
    SubSection,
)
from nomad_simulations.schema_packages.general import Simulation
from runschema.calculation import (
    ElectricFieldGradient as BaseElectricFieldGradient,
)
from runschema.calculation import (
    MagneticShielding as BaseMagneticShielding,
)
from runschema.calculation import (
    SpinSpinCoupling as BaseSpinSpinCoupling,
)

configuration = config.get_plugin_entry_point(
    'nomad_parser_magres.schema_packages:nomad_parser_magres_schema'
)

m_package = SchemaPackage()


class MagneticShielding(BaseMagneticShielding):
    """
    Section extensions for the Run.Calculation.MagneticShielding base section.
    """

    # ! These quantities should be implemented in `BaseMagneticShielding` as refs to the specific `AtomsState`

    m_def = Section(extends_base_section=True)

    atoms = Quantity(
        type=np.str_,
        shape=['n_atoms', 2],
        description="""
        Identifier for the atoms involved in the magnetic shielding tensor. This a list of
        `n_atoms` pairs of strings [atom_label, atom_index]. The atom index corresponds to the position
        on the list `System.atoms.labels`.
        """,
    )


class ElectricFieldGradient(BaseElectricFieldGradient):
    """
    Section extensions for the Run.Calculation.ElectricFieldGradient base section.
    """

    # ! These quantities should be implemented in `BaseElectricFieldGradient` as refs to the specific `AtomsState`

    m_def = Section(extends_base_section=True)

    atoms = Quantity(
        type=np.str_,
        shape=['n_atoms', 2],
        description="""
        Identifier for the atoms involved in the electric field gradient tensor. This a list of
        `n_atoms` pairs of strings [atom_label, atom_index]. The atom index corresponds to the position
        on the list `System.atoms.labels`.
        """,
    )


class SpinSpinCoupling(BaseSpinSpinCoupling):
    """
    Section extensions for the Run.Calculation.SpinspinCoupling base section.
    """

    # ! These quantities should be implemented in `BaseSpinSpinCoupling` as refs to the specific `AtomsState`g`

    m_def = Section(extends_base_section=True)

    atoms_1 = Quantity(
        type=np.str_,
        shape=['n_atoms', 2],
        description="""
        Identifier for the atoms involved in the spin-spin coupling J12 for the 1 atoms. This a list of
        `n_atoms` pairs of strings [atom_label, atom_index]. The atom index corresponds to the position
        on the list `System.atoms.labels`.
        """,
    )

    atoms_2 = Quantity(
        type=np.str_,
        shape=['n_atoms', 2],
        description="""
        Identifier for the atoms involved in the spin-spin coupling J12 for the 2 atoms. This a list of
        `n_atoms` pairs of strings [atom_label, atom_index]. The atom index corresponds to the position
        on the list `System.atoms.labels`.
        """,
    )


m_package.__init_metainfo__()
