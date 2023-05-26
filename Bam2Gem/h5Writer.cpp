/* Copyright (C) BGI-Reasearch - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by STOmics development team P_stomics_dev@genomics.cn, May 2017
 */
#include "h5Writer.h"

using namespace std;

static constexpr unsigned int version = 3;

H5Writer::H5Writer(string filename, unsigned int resolution)
    : maxexp(0), maxexon(0), min_x(INT_MAX), min_y(INT_MAX), max_x(0), max_y(0),
      resolution(resolution)
{
    m_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimsAttr[1] = { 1 };
    hid_t   attr;
    m_dataspace_id = H5Screate_simple(1, dimsAttr, NULL);
    attr     = H5Acreate(m_file_id, "version", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &version);
    H5Sclose(m_dataspace_id);
    H5Aclose(attr);

    m_group_id = H5Gcreate(m_file_id, "/geneExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_status   = H5Gclose(m_group_id);
    m_group_id =
        H5Gcreate(m_file_id, "/geneExp/bin1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    m_status = H5Gclose(m_group_id);
}

H5Writer::~H5Writer()
{
    m_status = H5Fclose(m_file_id);
}

bool H5Writer::store([[maybe_unused]] std::string& genename, [[maybe_unused]] exp_map& em)
{
    unsigned int coorX, coorY, geneid;

    // Incase the same genename in different contigs
    if (genes.count(genename) == 0)
    {
        genes[genename] = genePosList.size();
        // Store position data for current gene
        genePosList.emplace_back(genename.c_str(), dnbList.size(), em.size());

        // Store dnb data for current gene
        for (auto& [k, v] : em)
        {
            decodeKey(k, &coorX, &coorY, &geneid);
            dnbList.emplace_back(coorX, coorY, v.first);
            exonList.emplace_back(v.second);

            maxexp  = std::max(maxexp, v.first);
            maxexon = std::max(maxexon, v.second);
            min_x   = std::min(min_x, coorX);
            max_x   = std::max(max_x, coorX);
            min_y   = std::min(min_y, coorY);
            max_y   = std::max(max_y, coorY);
        }
    }
    else
    {
        auto  idx     = genes.at(genename);
        auto& genePos = genePosList[idx];

        vector< Dnb >          temp;
        vector< unsigned int > temp2;
        // Store dnb data for current gene
        for (auto& [k, v] : em)
        {
            decodeKey(k, &coorX, &coorY, &geneid);
            temp.emplace_back(coorX, coorY, v.first);
            temp2.emplace_back(v.second);

            maxexp  = std::max(maxexp, v.first);
            maxexon = std::max(maxexon, v.second);
            min_x   = std::min(min_x, coorX);
            max_x   = std::max(max_x, coorX);
            min_y   = std::min(min_y, coorY);
            max_y   = std::max(max_y, coorY);
        }

        // Insert the gene data with same name
        dnbList.insert(dnbList.begin() + genePos.offset + genePos.count, temp.begin(),
                       temp.end());
        exonList.insert(exonList.begin() + genePos.offset + genePos.count, temp2.begin(),
                        temp2.end());

        // Change the offset and count for genes behind the gene with same name
        genePos.count += em.size();
        for (unsigned int i = idx + 1; i < genePosList.size(); ++i)
            genePosList[i].offset += em.size();
    }

    return true;
}

bool H5Writer::dump()
{
    // Prepare datas
    std::transform(dnbList.begin(), dnbList.end(), dnbList.begin(), [&](Dnb& dnb) {
        Dnb res(dnb.x - min_x, dnb.y - min_y, dnb.cnt);
        return res;
    });

    // Create expression compound
    unsigned rank    = 1;
    hsize_t  dims[1] = { dnbList.size() };

    hid_t memtype, filetype;
    memtype  = H5Tcreate(H5T_COMPOUND, sizeof(Dnb));
    m_status = H5Tinsert(memtype, "x", HOFFSET(Dnb, Dnb::x), H5T_NATIVE_UINT);
    m_status = H5Tinsert(memtype, "y", HOFFSET(Dnb, Dnb::y), H5T_NATIVE_UINT);
    m_status = H5Tinsert(memtype, "count", HOFFSET(Dnb, Dnb::cnt), H5T_NATIVE_UINT);

    if (maxexp > USHRT_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 4);
        m_status = H5Tinsert(filetype, "x", 0, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "y", 4, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "count", 8, H5T_STD_U32LE);
    }
    else if (maxexp > UCHAR_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 2);
        m_status = H5Tinsert(filetype, "x", 0, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "y", 4, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "count", 8, H5T_STD_U16LE);
    }
    else
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 1);
        m_status = H5Tinsert(filetype, "x", 0, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "y", 4, H5T_STD_U32LE);
        m_status = H5Tinsert(filetype, "count", 8, H5T_STD_U8LE);
    }

    char expName[128] = { 0 };
    sprintf(expName, "/geneExp/bin1/expression");

    m_dataspace_id = H5Screate_simple(rank, dims, NULL);
    m_dataset_id   = H5Dcreate(m_file_id, expName, filetype, m_dataspace_id, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    m_status =
        H5Dwrite(m_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &*dnbList.begin());

    // Create expression attribute
    hsize_t dimsAttr[1] = { 1 };
    hid_t   attr;
    m_dataspace_id = H5Screate_simple(1, dimsAttr, NULL);
    attr     = H5Acreate(m_dataset_id, "minX", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &min_x);
    attr     = H5Acreate(m_dataset_id, "minY", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &min_y);
    attr     = H5Acreate(m_dataset_id, "maxX", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &max_x);
    attr     = H5Acreate(m_dataset_id, "maxY", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &max_y);
    attr = H5Acreate(m_dataset_id, "maxExp", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &maxexp);
    attr     = H5Acreate(m_dataset_id, "resolution", H5T_STD_U32LE, m_dataspace_id,
                     H5P_DEFAULT, H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &resolution);

    // Create expression of exonic
    char exonName[128] = { 0 };
    sprintf(exonName, "/geneExp/bin1/exon");

    m_dataspace_id = H5Screate_simple(rank, dims, NULL);
    if (maxexon > USHRT_MAX)
    {
        m_dataset_id = H5Dcreate(m_file_id, exonName, H5T_NATIVE_UINT, m_dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if (maxexon > UCHAR_MAX)
    {
        m_dataset_id = H5Dcreate(m_file_id, exonName, H5T_NATIVE_USHORT, m_dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        m_dataset_id = H5Dcreate(m_file_id, exonName, H5T_NATIVE_UCHAR, m_dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    m_status = H5Dwrite(m_dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        &*exonList.begin());

    m_dataspace_id = H5Screate_simple(1, dimsAttr, NULL);
    attr = H5Acreate(m_dataset_id, "maxExon", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
                     H5P_DEFAULT);
    m_status = H5Awrite(attr, H5T_NATIVE_UINT, &maxexon);

    // Create gene compound
    hid_t strtype;
    strtype  = H5Tcopy(H5T_C_S1);
    m_status = H5Tset_size(strtype, char_len);

    memtype  = H5Tcreate(H5T_COMPOUND, sizeof(GenePos));
    m_status = H5Tinsert(memtype, "gene", HOFFSET(GenePos, GenePos::gene), strtype);
    m_status =
        H5Tinsert(memtype, "offset", HOFFSET(GenePos, GenePos::offset), H5T_NATIVE_UINT);
    m_status =
        H5Tinsert(memtype, "count", HOFFSET(GenePos, GenePos::count), H5T_NATIVE_UINT);

    filetype = H5Tcreate(H5T_COMPOUND, char_len + 4 + 4);
    m_status = H5Tinsert(filetype, "gene", 0, strtype);
    m_status = H5Tinsert(filetype, "offset", char_len, H5T_STD_U32LE);
    m_status = H5Tinsert(filetype, "count", char_len + 4, H5T_STD_U32LE);

    dims[0]            = genePosList.size();
    char geneName[128] = { 0 };
    sprintf(geneName, "/geneExp/bin1/gene");
    m_dataspace_id = H5Screate_simple(rank, dims, NULL);
    m_dataset_id   = H5Dcreate(m_file_id, geneName, filetype, m_dataspace_id, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    m_status =
        H5Dwrite(m_dataset_id, filetype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &genePosList[0]);

    m_status = H5Tclose(strtype);
    m_status = H5Aclose(attr);
    m_status = H5Tclose(memtype);
    m_status = H5Tclose(filetype);
    m_status = H5Dclose(m_dataset_id);
    m_status = H5Sclose(m_dataspace_id);

    return true;
}
