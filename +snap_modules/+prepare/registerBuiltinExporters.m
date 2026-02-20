function exporters = registerBuiltinExporters()
% registerBuiltinExporters - Built-in channel exporters.

    exporters = snap_modules.plugins.prepare.bioformats_exporter_builtin();
end
