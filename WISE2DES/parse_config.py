from __future__ import print_function, division


def parse_config(config_file='wise2des.cfg', silent=False, debug=False):
    """

    read a config file

    might want to pass back the location of the config file

    """
    import os
    import ConfigParser
    config = ConfigParser.RawConfigParser()
    # read database connection info from config file
    # could check if existence in some default locations rather
    # using try

    if not silent:
        print('Reading config file', os.getcwd() + '/' + config_file)

    try:
        config.read(config_file)

    except Exception as e:
        print('Problem reading config file: ', config_file)
        print(e)
        print()
        print('config file', configfile, "does not exist")
        print('Now trying another location')
        config_file = os.path.join(os.environ["HOME"], '.config/', config_file)
        print('Trying ', configfile)
        cfg.read(config_file)
    except Exception as e:
        print('Problem reading config file: ', config_file)
        print('Now trying another location')
        config_file = os.path.join(os.environ["HOME"], '.config/', config_file)
        print('Trying ', configfile)
        cfg.read(config_file)

    # return db, host, user, password, table
    if debug:
        print()
        for section_name in config.sections():
            print('Section:', section_name)
            print('  Options:', config.options(section_name))
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
            print()
        print()

    return config


if __name__ == '__main__':

    config_file = 'wise2des.cfg'

    cfg = parse_config(config_file=config_file, debug=True)

    # help(cfg)
    tilename = cfg.get('des', 'tilename')

    print('TILENAME:', tilename)

    datapath = cfg.get('des', 'datapath')

    print('DATAPATH:', datapath)
