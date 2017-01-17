import os
import yaml
import warnings
import pkg_resources


def get_path_of_data_file(data_file):
    file_path = pkg_resources.resource_filename("gbm_drm_gen", 'data/%s' % data_file)

    return file_path


_config_file_name = 'gbmdrmgen_config.yml'


class Config(object):
    def __init__(self):

        # Read first the default configuration file
        default_configuration_path = get_path_of_data_file(_config_file_name)

        assert os.path.exists(default_configuration_path), \
            "Default configuration %s does not exist. Re-install gbm_drm_gen" % default_configuration_path

        with open(default_configuration_path) as f:

            try:

                self._configuration = yaml.safe_load(f)

            except:

                raise RuntimeError("Default configuration file %s cannot be parsed!" %
                                   (default_configuration_path))

            # This needs to be here for the _check_configuration to work




            self._default_path = default_configuration_path

        # Check if the user has a user-supplied config file under .threeML

        user_config_path = os.path.join(os.path.expanduser('~'), '.drm', _config_file_name)

        if os.path.exists(user_config_path):

            with open(user_config_path) as f:

                try:

                    self._configuration = yaml.safe_load(f)


                except:

                    raise RuntimeError("This file does not exist")

                self._filename = user_config_path

                print("Configuration read from %s" % (user_config_path))

        else:

            warnings.warn("Using default configuration from %s. "
                          "You might want to copy it to %s to customize it and avoid this warning."
                          % (self._default_path, user_config_path))

            # self._configuration = self._check_configuration(self._default_configuration_raw, self._default_path)
            self._filename = self._default_path

    def __getitem__(self, key):

        if key in self._configuration.keys():

            return self._configuration[key]

        else:

            raise ValueError("Configuration key %s does not exist in %s." % (key, self._filename))

    def __repr__(self):

        return yaml.dump(self._configuration, default_flow_style=False)


# Now read the config file, so it will be available as Config.c
gbm_drm_gen_config = Config()
