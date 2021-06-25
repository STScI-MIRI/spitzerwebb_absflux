import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.aperture import aperture_photometry
from astropy import wcs
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import simple_norm


def get_aper_flux(image, aperture, annulus_aperture, center_pix, chn_num):

    # IRAC correction coefficients
    a = np.zeros((6, 4))
    a[:, 0] = [1.0114, -3.536e-6, -6.826e-5, -1.618e-8, 1.215e-6, 1.049e-6]
    a[:, 1] = [1.0138, 8.401e-5, 3.345e-7, 1.885e-7, 1.438e-6, 1.337e-6]
    a[:, 2] = [1.0055, -3.870e-4, 4.600e-5, 1.956e-7, 2.078e-6, 9.970e-7]
    a[:, 3] = [1.0054, 2.332e-4, -8.234e-5, -1.881e-7, 6.520e-7, 9.415e-7]

    # IRAC correction coefficients (pixel phase)
    irac_phase_a = [0.0535, 0.0309]

    # source photometry
    phot_table = aperture_photometry(image, [aperture, annulus_aperture])

    # get a sigma clipped background measurement
    annulus_mask = annulus_aperture.to_mask(method="center")
    annulus_data = annulus_mask.multiply(image)
    mask = annulus_mask.data
    annulus_data_1d = annulus_data[mask > 0]
    med_bkg = np.nanmedian(annulus_data_1d)
    bkg = med_bkg * aperture.area

    ix_ref = center_pix[0]
    iy_ref = center_pix[1]
    corfac = (
        a[0, chn_num]
        + a[1, chn_num] * (ix_ref - 128.0)
        + a[2, chn_num] * (iy_ref - 128.0)
        + a[3, chn_num] * (ix_ref - 128.0) * (iy_ref - 128.0)
        + a[4, chn_num] * (ix_ref - 128.0) ** 2
        + a[5, chn_num] * (iy_ref - 128.0) ** 2
    )

    if (chn_num == 0) or (chn_num == 1):
        p_rad = np.sqrt(
            ((ix_ref - int(ix_ref)) - 0.5) + +2 + ((iy_ref - int(iy_ref)) - 0.5) ** 2
        )
        corfac2 = 1.0 + irac_phase_a[chn_num] * ((1.0 / np.sqrt(2.0 * np.pi)) - p_rad)
    else:
        corfac2 = 1.0

    # background subtracted flux
    return (phot_table["aperture_sum_0"][0] - bkg) * corfac * corfac2, corfac, corfac2


def disp_subimages(allfiles, sfluxes, aprads, new_center):
    # display the images
    n_rows = int(len(sfluxes) / 5) + 1
    fig = plt.figure(1, figsize=(12, 40))
    grid = ImageGrid(fig, 111, nrows_ncols=(n_rows, 5), axes_pad=0, share_all=False)

    i = 0
    for cfilename in allfiles:
        hdul = fits.open(cfilename)
        image = hdul[0].data
        w = wcs.WCS(hdul[0].header)

        center_pix = new_center.to_pixel(w)
        isize = image.shape
        if (0 < center_pix[0] < isize[0]) & (0 < center_pix[1] < isize[1]):

            mc = np.rint(center_pix).astype(int)

            subsize = 40
            image_cutout = Cutout2D(image, mc, [subsize, subsize], wcs=w)
            image = image_cutout.data
            center_pix = new_center.to_pixel(image_cutout.wcs)

            aperture = CircularAperture(center_pix, r=aprads[0])
            annulus_aperture = CircularAnnulus(
                center_pix, r_in=aprads[1], r_out=aprads[2]
            )

            norm = simple_norm(image, "sqrt", percent=99)
            grid[i].imshow(image, norm=norm, interpolation="nearest")

            aperture.plot(grid[i], color="white", lw=2, label="Photometry aperture")

            if sfluxes.mask[i]:
                dcolor = "red"
                ls = "dashed"
            else:
                dcolor = "magenta"
                ls = "dotted"

            annulus_aperture.plot(
                grid[i], color=dcolor, lw=2, linestyle=ls, label="Background annulus"
            )

            # grid[i].imshow(cimage, norm=norm_data, origin="lower")
            grid[i].axis("off")
            grid[i].set_xticks([])
            grid[i].set_yticks([])

            i += 1
