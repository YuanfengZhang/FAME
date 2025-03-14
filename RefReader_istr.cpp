//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#include <fstream>
#include <algorithm>

#include "RefReader_istr.h"

void readReference(const std::string& filename, std::vector<struct CpG>& cpgTab, std::vector<struct CpG>& cpgStartTab, std::vector<std::vector<char> >& genSeq, std::unordered_map<uint8_t, std::string>& chrMap, const bool humanOptFlag)
{
    std::string line;

    // Stores the sequence of each chromosome
    std::vector<char> seq;

    genSeq.reserve(MyConst::CHROMNUM);
    seq.reserve(MyConst::CHROMMAX);
    cpgTab.reserve(MyConst::CPGMAX);

    // Stores the chromosome index we are currently reading
    uint8_t chrIndex = 0;

    // Flag stating if we currently read a real chromosome assembly
    bool contFlag = false;
    // Flag stating if the last character read in the previous line was a 'c'
    bool lastC = false;

    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Opening genome reference file " << filename << " was unsuccessful! Exiting..." << std::endl;
        exit(1);
    }

    std::cout << "Start reading reference file " << filename << std::endl;

    uint32_t unIDCount = 1;

    while (getline(ifs, line))
    {
        // Test for id tag line
        if (line.begin() != line.end() && *(line.begin()) == '>' && (line.begin() + 1) != line.end())
        {
            // If we read primary assembly previously, write it to vectors
            if (contFlag)
            {
                seq.shrink_to_fit();
                genSeq.emplace_back(std::move(seq));
                // Reset buffer
                seq = std::vector<char>();
                seq.reserve(MyConst::CHROMMAX);
                // Reset flag for last 'c'
                lastC = false;
            }

            ++chrIndex;

            // Extract chromosome identifier from the header
            std::string chrHeader(line.begin() + 1, line.end());
            uint32_t pos = chrHeader.find_first_of(" \t");
            std::string chrID(chrHeader, 0, pos);

            // Check if the chromosome is valid
            if (humanOptFlag)
            {
                if (!isPrimaryHG(chrID))
                {
                    --chrIndex; // Decrement the chromosome index since this is not a valid chromosome
                    contFlag = false;
                    continue;
                }
            }

            // Ensure uniqueness of chromosome IDs
            for (const auto& chrP : chrMap)
            {
                if (chrP.second == chrID)
                {
                    std::cout << "WARNING: Chromosome identifier " << chrID << " found in header\n"
                              << chrHeader << "\nis not unique.";
                    chrID.append("_");
                    chrID.append(std::to_string(unIDCount++));
                    std::cout << "Renaming to " << chrID << "\n";
                    break;
                }
            }

            // Insert the chromosome ID into the map
            chrMap.insert(std::pair<uint8_t, std::string>(chrIndex - 1, chrID));
            contFlag = true;
            continue;
        }

        // Check if we are in a real primary chromosome assembly
        if (contFlag)
        {
            // Parse the line
            readLine(line, lastC, chrIndex, cpgTab, cpgStartTab, seq);
        }
        else
        {
            continue;
        }
    }

    // If we read primary assembly previously, write it to vectors
    if (contFlag)
    {
        seq.shrink_to_fit();
        genSeq.emplace_back(std::move(seq));
    }

    cpgTab.shrink_to_fit();
    genSeq.shrink_to_fit();
    ifs.close();

    std::cout << "Done reading reference file" << std::endl;
}

bool isPrimaryHG(std::string chrID)
{
    // Include all standard chromosomes (chr1 to chr22)
    for (unsigned int i = 1; i <= 22; ++i)
    {
        if (chrID == ("chr" + std::to_string(i)))
        {
            return true;
        }
    }

    // Include additional special chromosomes
    const std::vector<std::string> specialChromosomes = {"chrX", "chrY", "chrM", "chrMT", "lambda", "pUC19"};
    if (std::find(specialChromosomes.begin(), specialChromosomes.end(), chrID) != specialChromosomes.end())
    {
        return true;
    }

    return false;
}