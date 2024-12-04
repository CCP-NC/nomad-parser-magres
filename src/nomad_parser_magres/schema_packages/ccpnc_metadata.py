import numpy as np
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import JSON, Quantity, SubSection


class MaterialProperties(ArchiveSection):
    # Note from @JosePizarro3: note we have all these information somewhere else in the `nomad_simulations` schema.
    # Nevertheless, if you feel it is better to keep these quantities here for clarity, it is totally fine for me.
    chemical_name = Quantity(
        type=str,
        description="""
        Free-text chemical name assigned by users.
        """,
    )

    chemical_name_tokens = Quantity(
        type=str,
        shape=["*"],
        description="""
        Free-text chemical name, but tokenised to take individual words in the name to assist in wildcard searches.
        """,
    )

    formula = Quantity(
        # type=[(str, int)],  # better to use JSON type
        type=JSON,
        shape=["*"],
        description="""
        Dictionary containing the species (chemical symbol of an element in the material) as keys and
        number of atoms of that element in the material as their value.
        """,
    )

    stoichiometry = Quantity(
        type=JSON,
        shape=["*"],
        description="""
        Reduced proportion of materials details.
        """,
    )

    elements_ratios = Quantity(
        type=np.float64,
        shape=["*"],
        description="""
        Ratio of constituent elements (each element is a number between 0 and 1).
        """,
    )

    # Note from @JosePizarro3: in the `nomad_simulations` schema we have a sub-section under `archive.data.model_system[*].chemical_formula`
    # where we compiled a bunch of different formats. There, `chemical_formula.descriptive` is selected depending on the
    # specific case (in organic and inorganic chemistry, the descriptive formula is different).
    chemical_formula_descriptive = Quantity(
        type=str,
        description="""
        Formula as a string, e.g., 'C2H6O'.
        """,
    )


class ORCID(ArchiveSection):
    orcid_id = Quantity(
        type=JSON,
        shape=["*"],
        description="""
        Dictionary containing the ORCID IDs of the author (keys) and uploader (values) profiles.
        """,
    )


class CCPNCRecord(ArchiveSection):
    visible = Quantity(
        type=bool,
        description="""
        A boolean value that indicates if the record is to be hidden or available to be returned when searched
        """,
    )

    immutable_id = Quantity(
        type=str,
        description="""
        7 digit unique record identifier.
        """,
    )


class ExternalDatabaseReference(ArchiveSection):
    external_database_name = Quantity(
        type=str,
        description="""
        External database name where additional information on the material exists
        """,
    )

    external_database_reference_code = Quantity(
        type=str,
        description="""
        Specific database code pointing to the material or a polymorphic form of the material.
        """,
    )


class FreeTextMetadata(ArchiveSection):
    uploader_author_notes = Quantity(
        type=str,
        description="""
        Additional metadata that authors want to indicate about the computation.
        """,
    )

    structural_descriptor_notes = Quantity(
        type=str,
        description="""
        Additional notes specific to the polymorphic forms of the material.
        """,
    )


class CCPNCMetadata(ArchiveSection):
    material_properties = SubSection(section_def=MaterialProperties)
    orcid = SubSection(section_def=ORCID)
    ccpnc_record = SubSection(section_def=CCPNCRecord)
    external_database_reference = SubSection(section_def=ExternalDatabaseReference)
    free_text_metadata = SubSection(section_def=FreeTextMetadata)
