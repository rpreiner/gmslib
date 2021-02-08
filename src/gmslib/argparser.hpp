//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2019
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------


#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>


class ArgParser
{
private:
	struct Arg
	{
		std::string name;
		std::string desc;
		std::string value;
	};

	std::vector<Arg> arglist;
	std::string noteStr = "";

public:
	bool matchArgs(int argc, char** argv)
	{
		for (int i = 1; i < argc - 1; i++)
		{
			bool found = false;
			for (auto& arg : arglist)
				if (argv[i] == "-" + arg.name)
				{
					arg.value = argv[i + 1];
					i++;
					found = true;
					break;
				}

			// Argument not found
			if (!found)
			{
				std::cerr << "Unknown argument '" << argv[i] << "'" << std::endl;
				return false;
			}
		}
		return true;
	}

	void addNote(const std::string& note)
	{
		noteStr = note;
	}

	void addArgument(const std::string& name, const std::string& description)
	{
		Arg arg;
		arg.name = name;
		arg.desc = description;
		arg.value = "";
		
		arglist.push_back(arg);
	}

	const std::string getArgument(const std::string& name) const
	{
		for (auto& arg : arglist)
			if (arg.name == name)
				return arg.value;
		
		return "";
	}

	bool getBool(const std::string& name, bool defaultValue) const
	{
		for (auto& arg : arglist)
			if (arg.name == name && arg.value != "")
			{
				std::string s;
				for (char c : arg.value) s += tolower(c);
				
				if (s == "true" || arg.value == "1")
					return true;
				else if (s == "false" || arg.value == "0")
					return false;
				else
				{
					std::cerr << "Invalid bool argument '" << arg.value << "'" << std::endl;
					exit(1);
				}
			}
				
		return defaultValue;
	}

	float getFloat(const std::string& name, float defaultValue) const
	{
		for (auto& arg : arglist)
			if (arg.name == name && arg.value != "")
			{
				try {
					float value = stof(arg.value);
					return value;
				}
				catch (std::invalid_argument e) {
					std::cerr << "Invalid float argument '" << arg.value << "'" << std::endl;
					exit(1); 
				}
			}

		return defaultValue;
	}

	unsigned getUint (const std::string& name, unsigned defaultValue) const
	{
		for (auto& arg : arglist)
			if (arg.name == name && arg.value != "")
			{
				try {
					unsigned value = stoi(arg.value);
					return value;
				}
				catch (std::invalid_argument e) {
					std::cerr << "Invalid unsigned int argument '" << arg.value << "'" << std::endl;
					exit(1);
				}
			}

		return defaultValue;
	}
	
	std::string helpStr() const
	{
		std::stringstream s;
		s << std::left << "\n";
		s << std::setw(14) << "Switch"
		  << std::setw(70) << "Description"
		  << "\n\n";
		
		for (auto& arg : arglist)
			s << std::setw(14) << "-"+arg.name << std::setw(70) << arg.desc << "\n";
		s << "\n";

		if (noteStr != "")
			s << noteStr << "\n";

		return s.str();
	}
};
