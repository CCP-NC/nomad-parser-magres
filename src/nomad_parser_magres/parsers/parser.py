import json
import os
from typing import TYPE_CHECKING, Optional

import numpy as np

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import EntryArchive
    from nomad_simulations.schema_packages.model_system import Cell
    from structlog.stdlib import BoundLogger

from typing import Dict

from nomad.app.v1.models.models import MetadataRequired
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Link, TaskReference
from nomad.parsing import MatchingParser
from nomad.parsing.file_parser import Quantity, TextParser
from nomad.search import search
from nomad.units import ureg
from nomad.utils import extract_section
from nomad_simulations.schema_packages.atoms_state import AtomsState
from nomad_simulations.schema_packages.general import Program
from nomad_simulations.schema_packages.model_method import (
    DFT,
    ModelMethod,
    XCFunctional,
)
from nomad_simulations.schema_packages.model_system import AtomicCell, ModelSystem
from nomad_simulations.schema_packages.numerical_settings import KMesh, KSpace

# utility function used to get auxiliary files next to the `mainfile`
from nomad_parser_magres.parsers.utils import get_files
from nomad_parser_magres.schema_packages.ccpnc_metadata import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
)
from nomad_parser_magres.schema_packages.package import CCPNCSimulation as Simulation
from nomad_parser_magres.schema_packages.package import (
    ElectricFieldGradient,
    ElectricFieldGradients,
    MagneticShieldingTensor,
    MagneticSusceptibility,
    Outputs,
    SpinSpinCoupling,
)
from nomad_parser_magres.schema_packages.workflow import (
    NMRMagRes,
    NMRMagResMethod,
    NMRMagResResults,
)

re_float = r" *[-+]?\d+\.\d*(?:[Ee][-+]\d+)? *"

configuration = config.get_plugin_entry_point(
    "nomad_parser_magres.parsers:nomad_parser_magres_plugin"
)


class MagresFileParser(TextParser):
    def __init__(self):
        super().__init__()

    def init_quantities(self):
        self._quantities = [
            Quantity("lattice_units", r"units *lattice *([a-zA-Z]+)"),
            Quantity("atom_units", r"units *atom *([a-zA-Z]+)"),
            Quantity("ms_units", r"units *ms *([a-zA-Z]+)"),
            Quantity("efg_units", r"units *efg *([a-zA-Z]+)"),
            Quantity("efg_local_units", r"units *efg_local *([a-zA-Z]+)"),
            Quantity("efg_nonlocal_units", r"units *efg_nonlocal *([a-zA-Z]+)"),
            Quantity("isc_units", r"units *isc *([a-zA-Z\^\d\.\-]+)"),
            Quantity("isc_fc_units", r"units *isc_fc *([a-zA-Z\^\d\.\-]+)"),
            Quantity("isc_spin_units", r"units *isc_spin *([a-zA-Z\^\d\.\-]+)"),
            Quantity(
                "isc_orbital_p_units", r"units *isc_orbital_p *([a-zA-Z\^\d\.\-]+)"
            ),
            Quantity(
                "isc_orbital_d_units", r"units *isc_orbital_d *([a-zA-Z\^\d\.\-]+)"
            ),
            Quantity("sus_units", r"units *sus *([a-zA-Z\^\d\.\-]+)"),
            Quantity("cutoffenergy_units", r"units *calc\_cutoffenergy *([a-zA-Z]+)"),
            Quantity(
                "calculation",
                r"([\[\<]*calculation[\>\]]*[\s\S]+?)(?:[\[\<]*\/calculation[\>\]]*)",
                sub_parser=TextParser(
                    quantities=[
                        Quantity("code", r"calc\_code *([a-zA-Z]+)"),
                        Quantity(
                            "code_version", r"calc\_code\_version *([a-zA-Z\d\.]+)"
                        ),
                        Quantity(
                            "code_hgversion",
                            r"calc\_code\_hgversion ([a-zA-Z\d\:\+\s]*)\n",
                            flatten=False,
                        ),
                        Quantity(
                            "code_platform", r"calc\_code\_platform *([a-zA-Z\d\_]+)"
                        ),
                        Quantity("name", r"calc\_name *([\w]+)"),
                        Quantity("comment", r"calc\_comment *([\w]+)"),
                        Quantity("xcfunctional", r"calc\_xcfunctional *([\w]+)"),
                        Quantity(
                            "cutoffenergy",
                            rf"calc\_cutoffenergy({re_float})(?P<__unit>\w+)",
                        ),
                        Quantity(
                            "pspot",
                            r"calc\_pspot *([\w]+) *([\w\.\|\(\)\=\:]+)",
                            repeats=True,
                        ),
                        Quantity(
                            "kpoint_mp_grid",
                            r"calc\_kpoint\_mp\_grid *([\w]+) *([\w]+) *([\w]+)",
                        ),
                        Quantity(
                            "kpoint_mp_offset",
                            rf"calc\_kpoint\_mp\_offset({re_float * 3})$",
                        ),
                    ]
                ),
            ),
            Quantity(
                "atoms",
                r"([\[\<]*atoms[\>\]]*[\s\S]+?)(?:[\[\<]*\/atoms[\>\]]*)",
                sub_parser=TextParser(
                    quantities=[
                        Quantity("lattice", rf"lattice({re_float * 9})"),
                        Quantity("symmetry", r"symmetry *([\w\-\+\,]+)", repeats=True),
                        Quantity(
                            "atom",
                            rf"atom *([a-zA-Z]+) *[a-zA-Z\d]* *([\d]+) *({re_float * 3})",
                            repeats=True,
                        ),
                    ]
                ),
            ),
            Quantity(
                "magres",
                r"([\[\<]*magres[\>\]]*[\s\S]+?)(?:[\[\<]*\/magres[\>\]]*)",
                sub_parser=TextParser(
                    quantities=[
                        Quantity(
                            "ms", rf"ms *(\w+) *(\d+)({re_float * 9})", repeats=True
                        ),
                        Quantity(
                            "efg", rf"efg *(\w+) *(\d+)({re_float * 9})", repeats=True
                        ),
                        Quantity(
                            "efg_local",
                            rf"efg_local *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "efg_nonlocal",
                            rf"efg_nonlocal *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "isc",
                            rf"isc *(\w+) *(\d+) *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "isc_fc",
                            rf"isc_fc *(\w+) *(\d+) *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "isc_orbital_p",
                            rf"isc_orbital_p *(\w+) *(\d+) *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "isc_orbital_d",
                            rf"isc_orbital_d *(\w+) *(\d+) *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity(
                            "isc_spin",
                            rf"isc_spin *(\w+) *(\d+) *(\w+) *(\d+)({re_float * 9})",
                            repeats=True,
                        ),
                        Quantity("sus", rf"sus *({re_float * 9})", repeats=True),
                    ]
                ),
            ),
        ]


class CcpncMagresParser(MagresParser):


    def parse_json_file(self, filepath: str, logger: "BoundLogger") -> Optional[CCPNCMetadata]:
        """Parse the JSON file and extract relevant information."""
        magres_json_file = get_files(
            pattern="MRD*.json", filepath=filepath, stripname=self.basename
        )
        if not magres_json_file:
            logger.warning("No JSON file found.")
            return None
        with open(magres_json_file[0]) as f:
            magres_json_data = json.load(f)
        ccpnc_metadata = CCPNCMetadata()
        material_properties = MaterialProperties()
        orcid = ORCID()
        ccpnc_record = CCPNCRecord()
        external_database_reference = ExternalDatabaseReference()
        free_text_metadata = FreeTextMetadata()

        material_properties.chemical_name = magres_json_data.get("chemname", "")
        material_properties.formula = magres_json_data.get("formula", "")
        material_properties.stoichiometry = magres_json_data.get("stochiometry", "")
        material_properties.elements_ratios = magres_json_data.get("elements_ratios", "")
        # material_properties.chemical_name_tokens =
        orcid.orcid_id = magres_json_data.get("ORCID", "")
        # ccpnc_record.visible =
        ccpnc_record.immutable_id = magres_json_data.get("immutable_id", "")
        version_metadata = magres_json_data.get("version_metadata", {})
        external_database_reference.external_database_name = version_metadata.get("extref_type", "")
        external_database_reference.external_database_reference_code = version_metadata.get("extref_code", "")
        free_text_metadata.uploader_author_notes = version_metadata.get("notes", "")
        free_text_metadata.structural_descriptor_notes = version_metadata.get("chemform", "")

        ccpnc_metadata.material_properties = material_properties
        ccpnc_metadata.orcid = orcid
        ccpnc_metadata.ccpnc_record = ccpnc_record
        ccpnc_metadata.external_database_reference = external_database_reference
        ccpnc_metadata.free_text_metadata = free_text_metadata
        return ccpnc_metadata

    def parse(
        self,
        filepath: str,
        archive: "EntryArchive",
        logger: "BoundLogger",
        child_archives: Dict[str, EntryArchive] = None,
    ) -> None:


        # Populate `CCPNCMetadata` (note the `pattern` has to match the aux file generated by the MongoDB CCP-NC)
        magres_json_file = get_files(
            pattern="magres*.json", filepath=self.mainfile, stripname=self.basename
        )
        if magres_json_file is not None:
            ccpnc_metadata = CCPNCMetadata()
            # TODO: populate `ccpnc_metadata` model from `magres_json_file` HERE
            # ...
            # ...
            simulation.ccpnc_metadata = ccpnc_metadata
