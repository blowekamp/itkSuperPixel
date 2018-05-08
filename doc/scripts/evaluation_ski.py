#!/usr/bin/env python

import SimpleITK as sitk
import skimage as ski
import skimage.segmentation
import numpy as np
import timeit


def mask_label_contour(image, seg):
   """Combine an image and segmentation by masking the segmentation contour.

   For an input image (scalar or vector), and a multi-label
   segmentation image, creates an output image where the countour of
   each label masks the input image to black."""
   return sitk.Mask(image, sitk.LabelContour(seg+1)==0)

# this script generates images to compare ski-image SLIC
# implementaiton vs ours.

# We have slightly different parameterizations. The image is 512x512,
# if we target 256 superpixels of size 32x32 we have simular
# parameters for each implementation.

img=sitk.ReadImage("/home/blowekamp/src/scikit-image/skimage/data/astronaut.png")

aimg_lab=ski.color.rgb2lab(sitk.GetArrayFromImage(img))
ski_slic_aimg=skimage.segmentation.slic(aimg_lab,n_segments=256,convert2lab=False)

sitk.WriteImage(mask_label_contour(img, sitk.GetImageFromArray(ski_slic_aimg))
, "astronaut_ski_slic.png")

print(min(timeit.repeat(lambda: skimage.segmentation.slic(aimg_lab,n_segments=256,convert2lab=False), number=1, repeat=5)))

img_lab = sitk.GetImageFromArray(aimg_lab, isVector=True)
sitk_slic_img=sitk.SLIC(img_lab, [32,32], maximumNumberOfIterations=10)
sitk.WriteImage(mask_label_contour(img, sitk_slic_img), "astronaut_sitk_slic.png")

print(min(timeit.repeat(lambda: sitk.SLIC(img_lab, [32,32], maximumNumberOfIterations=10), number=1, repeat=5)))







