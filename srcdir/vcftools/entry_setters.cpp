/*
	entry_setters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "entry.h"

void variant_file_entry::set_CHROM(const string &in)
{
	CHROM = in;
}

void variant_file_entry::set_POS(const int in)
{
	POS = in;
}

void variant_file_entry::set_ID(const string &in)
{
	ID = in;
}

void variant_file_entry::set_REF(const string &in)
{
	REF = in;
}

void variant_file_entry::add_ALT_allele(const string &in)
{
	if (in != ".")
	{
		if (find(ALT.begin(), ALT.end(),in) == ALT.end())
		{
			ALT.push_back(in);
		}
		else
		{
			LOG.warning(" Duplicate alternate alleles found at " + CHROM + ":" + LOG.int2str(POS));
			ALT.push_back(in);
		}
	}
	parsed_ALT = true;
}

void variant_file_entry::add_FILTER_entry(const string &in)
{
	if ((in != "PASS") && (in != "."))
		if (find(FILTER.begin(), FILTER.end(), in) == FILTER.end())
			FILTER.push_back(in);
	sort(FILTER.begin(), FILTER.end());
	parsed_FILTER = true;
}
