import numpy as np
from nomad.metainfo import Quantity, Section, SubSection
from nomad.datamodel.data import ArchiveSection

class MaterialProperties(ArchiveSection):
    chemical_name = Quantity(
        type=str,
        shape=[],
        description="Free-text chemical name assigned by users"
    )

    chemical_name_tokens = Quantity(
        type=str,
        shape=['*'],
        description="Free-text chemical name, but tokenised to take individual words in the name to assist in wildcard searches"
    )

    formula = Quantity(
        type=[(str, int)],
        shape=['*'],
        description="(species, n) array where species refers to the chemical symbol of an element in the material, n denotes the number of atoms of that element in the material"
    )

    stoichiometry = Quantity(
        type=[(str, int)],
        shape=['*'],
        description="Reduced proportion of materials details"
    )

    elements_ratios = Quantity(
        type=np.dtype(np.float64),
        shape=['*'],
        description="Ratio of constituent elements (each element is a number between 0 and 1)"
    )

    chemical_formula_descriptive = Quantity(
        type=str,
        shape=[],
        description="Formula as a string eg. C2H6O"
    )

class ORCID(ArchiveSection):
    orcid_id = Quantity(
        type=[(str, str)],
        shape=['*'],
        description="ORCID ID of author and uploader profiles"
    )

class CCPNCRecord(ArchiveSection):
    visible = Quantity(
        type=bool,
        shape=[],
        description="A boolean value that indicates if the record is to be hidden or available to be returned when searched"
    )

    immutable_id = Quantity(
        type=str,
        shape=[],
        description="7 digit unique record identifier"
    )

class ExternalDatabaseReference(ArchiveSection):
    external_database_name = Quantity(
        type=str,
        shape=[],
        description="External database name where additional information on the material exists"
    )

    external_database_reference_code = Quantity(
        type=str,
        shape=[],
        description="Specific database code pointing to the material or a polymorphic form of the material"
    )

class FreeTextMetadata(ArchiveSection):
    uploader_author_notes = Quantity(
        type=str,
        shape=[],
        description="Additional metadata that authors want to indicate about the computation"
    )

    structural_descriptor_notes = Quantity(
        type=str,
        shape=[],
        description="Additional notes specific to the polymorphic forms of the material"
    )

class CCPNCMetadata(ArchiveSection):
    material_properties = SubSection(section_def=MaterialProperties)
    orcid = SubSection(section_def=ORCID)
    ccpnc_record = SubSection(section_def=CCPNCRecord)
    external_database_reference = SubSection(section_def=ExternalDatabaseReference)
    free_text_metadata = SubSection(section_def=FreeTextMetadata)