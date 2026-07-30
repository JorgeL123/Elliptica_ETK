#include <vector>
#include <string>

#include <elliptica_id_reader.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

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

        vel[i]               = idr->field[i_vx][i];
        vel[i + npoints]     = idr->field[i_vy][i];
        vel[i + 2*npoints]   = idr->field[i_vz][i];
    }

    elliptica_id_reader_free(idr);

    CCTK_INFO("Elliptica initial data successfully imported");
}
