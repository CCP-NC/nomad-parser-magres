from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class MagresSchemaPackageEntryPoint(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_parser_magres.schema_packages.package import m_package

        return m_package


nomad_parser_magres_schema = MagresSchemaPackageEntryPoint(
    name='MagresSchemaPackageEntryPoint',
    description='Entry point for the Magres code-specific schema.',
)
