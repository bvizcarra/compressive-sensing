#ifndef __CONFIG_FILE_H__
#define __CONFIG_FILE_H__

#include <string>
#include <map>

#include "Chameleon.h"

/// ConfigFile is a class for reading config files using the class Chamelon.
/// It was written by Rene Nyffenegger.
/// It is available at the website http://www.adp-gmbh.ch/cpp/config_file.html.
class ConfigFile {
	std::map<std::string,Chameleon> content_;

public:
	ConfigFile(std::string const& configFile);

	Chameleon const& Value(std::string const& section, std::string const& entry) const;

	Chameleon const& Value(std::string const& section, std::string const& entry, double value);
	Chameleon const& Value(std::string const& section, std::string const& entry, std::string const& value);
};

#endif