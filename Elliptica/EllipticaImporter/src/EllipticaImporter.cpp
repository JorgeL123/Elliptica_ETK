#include <vector>
#include <string>

#include <elliptica_id_reader_lib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C"
void EllipticaImporter_check_parameters(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_PARAMETERS;

    if (CCTK_EQUALS(checkpoint_path, ""))
    {
        CCTK_PARAMWARN("EllipticaImporter::checkpoint_path must be set to "
                        "the full path of an Elliptica checkpoint file");
    }

    if (!(CCTK_EQUALS(type, "BBH")  ||
          CCTK_EQUALS(type, "BH")   ||
          CCTK_EQUALS(type, "BNS")  ||
          CCTK_EQUALS(type, "NS")   ||
          CCTK_EQUALS(type, "BHNS")))
    {
        CCTK_PARAMWARN("EllipticaImporter::type must be one of "
                        "\"BBH\", \"BH\", \"BNS\", \"NS\", \"BHNS\"");
    }
}

extern "C"
void EllipticaImporter(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    CCTK_INFO("Setting up Elliptica initial data");

    const int npoints =
        cctk_lsh[0] *
        cctk_lsh[1] *
        cctk_lsh[2];

    // Stride between vector components (e.g. vel[3]) in Cactus's flat
    // storage is fixed by the ALLOCATED array shape, not the logical
    // local shape -- these can differ (e.g. SIMD/vector padding), so
    // this must NOT be computed from cctk_lsh.
    const int np_alloc =
        cctk_ash[0] *
        cctk_ash[1] *
        cctk_ash[2];

    // Coordinate arrays expected by Elliptica
    std::vector<CCTK_REAL> xx(npoints);
    std::vector<CCTK_REAL> yy(npoints);
    std::vector<CCTK_REAL> zz(npoints);

    #pragma omp parallel for
    for (int i = 0; i < npoints; ++i)
    {
        xx[i] = x[i];
        yy[i] = y[i];
        zz[i] = z[i];
    }

    CCTK_INFO("Initializing Elliptica ID Reader");

    if (CCTK_EQUALS(checkpoint_path, ""))
    {
        CCTK_ERROR("EllipticaImporter::checkpoint_path is not set");
    }

    Elliptica_ID_Reader_T *idr =
        elliptica_id_reader_init(
            checkpoint_path,
            "generic_MT_safe"
        );

    if (idr == nullptr)
    {
        CCTK_ERROR("Could not initialize Elliptica ID Reader");
    }

    idr->ifields =
        "alpha,betax,betay,betaz,"
        "adm_gxx,adm_gxy,adm_gxz,"
        "adm_gyy,adm_gyz,adm_gzz,"
        "adm_Kxx,adm_Kxy,adm_Kxz,"
        "adm_Kyy,adm_Kyz,adm_Kzz,"
        "grhd_rho,grhd_epsl,grhd_p,"
        "grhd_vx,grhd_vy,grhd_vz";

    idr->npoints = npoints;

    idr->x_coords = xx.data();
    idr->y_coords = yy.data();
    idr->z_coords = zz.data();

    // Recommended settings for BHNS/BNS
    idr->set_param(
        "BH_filler_method",
        "ChebTn_Ylm_perfect_s2",
        idr
    );

    idr->set_param(
        "ADM_B1I_form",
        "zero",
        idr
    );

    CCTK_INFO("Interpolating Elliptica initial data");

    elliptica_id_reader_interpolate(idr);

    // Cache field indices
    const int i_alpha = idr->indx("alpha");
    const int i_betax = idr->indx("betax");
    const int i_betay = idr->indx("betay");
    const int i_betaz = idr->indx("betaz");

    const int i_gxx = idr->indx("adm_gxx");
    const int i_gxy = idr->indx("adm_gxy");
    const int i_gxz = idr->indx("adm_gxz");
    const int i_gyy = idr->indx("adm_gyy");
    const int i_gyz = idr->indx("adm_gyz");
    const int i_gzz = idr->indx("adm_gzz");

    const int i_Kxx = idr->indx("adm_Kxx");
    const int i_Kxy = idr->indx("adm_Kxy");
    const int i_Kxz = idr->indx("adm_Kxz");
    const int i_Kyy = idr->indx("adm_Kyy");
    const int i_Kyz = idr->indx("adm_Kyz");
    const int i_Kzz = idr->indx("adm_Kzz");

    const int i_rho   = idr->indx("grhd_rho");
    const int i_eps   = idr->indx("grhd_epsl");
    const int i_press = idr->indx("grhd_p");

    const int i_vx = idr->indx("grhd_vx");
    const int i_vy = idr->indx("grhd_vy");
    const int i_vz = idr->indx("grhd_vz");

    const int indices[] = {
        i_alpha, i_betax, i_betay, i_betaz,
        i_gxx, i_gxy, i_gxz, i_gyy, i_gyz, i_gzz,
        i_Kxx, i_Kxy, i_Kxz, i_Kyy, i_Kyz, i_Kzz,
        i_rho, i_eps, i_press, i_vx, i_vy, i_vz
    };
    for (int k = 0; k < 22; ++k)
    {
        if (indices[k] < 0)
        {
            CCTK_ERROR("EllipticaImporter: idr->indx() returned -1 for one "
                        "or more requested fields -- checkpoint may not "
                        "contain all fields listed in idr->ifields for "
                        "this run's type. Aborting before writing to "
                        "grid functions.");
        }
    }

    CCTK_INFO("Copying Elliptica data to Cactus grid functions");

    #pragma omp parallel for
    for (int i = 0; i < npoints; ++i)
    {
        // ADM variables
        alp[i]   = idr->field[i_alpha][i];

        betax[i] = idr->field[i_betax][i];
        betay[i] = idr->field[i_betay][i];
        betaz[i] = idr->field[i_betaz][i];

        gxx[i] = idr->field[i_gxx][i];
        gxy[i] = idr->field[i_gxy][i];
        gxz[i] = idr->field[i_gxz][i];
        gyy[i] = idr->field[i_gyy][i];
        gyz[i] = idr->field[i_gyz][i];
        gzz[i] = idr->field[i_gzz][i];

        kxx[i] = idr->field[i_Kxx][i];
        kxy[i] = idr->field[i_Kxy][i];
        kxz[i] = idr->field[i_Kxz][i];
        kyy[i] = idr->field[i_Kyy][i];
        kyz[i] = idr->field[i_Kyz][i];
        kzz[i] = idr->field[i_Kzz][i];

        // Hydrodynamic variables
        rho[i]   = idr->field[i_rho][i];
        eps[i]   = idr->field[i_eps][i];
        press[i] = idr->field[i_press][i];

        // HydroBase::vel is declared as `CCTK_REAL vel[3] type = GF`: one
        // flat buffer of 3 * np_alloc elements, components separated by
        // the ALLOCATED array size (np_alloc), not the logical size
        // (npoints). This matches the convention used in HydroBase's own
        // Initialization.c (HydroBase_Zero).
        vel[i]               = idr->field[i_vx][i];
        vel[i + np_alloc]    = idr->field[i_vy][i];
        vel[i + 2*np_alloc]  = idr->field[i_vz][i];
    }

    elliptica_id_reader_free(idr);

    CCTK_INFO("Elliptica initial data successfully imported");
}

    CCTK_INFO("Elliptica initial data successfully imported");
}
