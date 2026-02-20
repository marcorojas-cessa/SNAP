function providers = registerBuiltinProviders()
% registerBuiltinProviders - Built-in image-library reader providers.

    providers = snap_modules.plugins.prepare.bioformats_reader_builtin();
end
