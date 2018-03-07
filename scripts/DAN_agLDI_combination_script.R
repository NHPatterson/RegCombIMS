#Code for the use of the IMS data analysis pipeline described in:
#"Co-registration and analysis of multiple imag-ing mass spectrometry datasets targeting dif-ferent analytes"
#Nathan Heath Patterson, Ethan Yang, Elizabeth-Ann Kranjec, Pierre Chaurand
#...

#install dependencies listed below:

require(Cardinal) #see cardinalmsi.org for instruction on installation
require(RegCombIMS)
require(jpeg)
require(RNiftyReg)
require(abind)
require(plyr)

#load:

#get IMS template images
d_reg_img_pca <- resultToArray(d001_cer, d_pca, idx=1)
s_reg_img_pca <- resultToArray(s000_cer, s_pca, idx=1)

##import image for registration, (he2 is loaded in the data namespace but code below shows how to load image, greyscale it and pad it for registration wtih IMS):
# he = readJPEG('./H_E_stain/crop_cerebellum_resize.jpg')
# 
# ##image file to greyscale for registration & transpose to match IMS data:
# he <- apply(he,1:2,mean)
# he <- t(he)
# 
# ##clean up borders of image..
# he[he > 0.95] = 0
# 
# ##pad borders of image with zeroes, if IMS data is registered and some is outside the coordinates of the H&E image, it will be cut-off
# he2 <-matrix(0,ncol=ncol(he) + 60,nrow= nrow(he) + 60)
# he2[30:(nrow(he)+29), 30:(ncol(he)+29)] <- he

#do affine registration, initialization
reg_d_to_he_affine_pca <- niftyreg(d_reg_img_pca, he2, scope=c("affine"),verbose=FALSE)
reg_s_to_he_affine_pca <- niftyreg(s_reg_img_pca, he2, scope=c("affine"),verbose=FALSE)

image(he2, main= "H&E")
image(reg_d_to_he_affine_pca$image, main="Affine_DAN_pca")
image(he2, main= "H&E")
image(reg_s_to_he_affine_pca$image, main= "Affine_AgLDI_pca")

#do non-linear registration
reg_d_to_he_nonlinear_pca <- niftyreg(d_reg_img_pca, he2, scope=c("nonlinear"),init=forward(reg_d_to_he_affine_pca), maxIterations = 1000,nBins = 128L)
reg_s_to_he_nonlinear_pca <- niftyreg(s_reg_img_pca, he2, scope=c("nonlinear"),init=forward(reg_s_to_he_affine_pca), maxIterations = 1000,nBins = 128L)

image(he2, main= "H&E")
image(reg_d_to_he_nonlinear_pca$image, main="Non-linear_DAN_pca")
image(he2, main= "H&E")
image(reg_s_to_he_nonlinear_pca$image, main= "Non-linear_AgLDI_pca")

##combine datasets: 

#transform the datasets using their respective transformations
tformed_d001_cer = applyNifty(d001_cer, reg_d_to_he_nonlinear_pca, padding=0, interpolation = 0)
tformed_s000_cer = applyNifty(s000_cer, reg_s_to_he_nonlinear_pca, padding=0, interpolation = 0)

image(d001_cer, mz=885.54, main="(-)DAN pre-transform")
image(tformed_d001_cer, mz=885.54, main="(-)DAN non-linear transformed")
image(s000_cer, mz=493.25, main="AgLDI pre-transform")
image(tformed_s000_cer, mz=493.25, main="AgLDI non-linear transformed")

##combine the datasets into a single dataset:
#rescale intensity data between 0 and 1 for each dataset
iData(tformed_d001_cer) = iData(tformed_d001_cer) / max(iData(tformed_d001_cer))
iData(tformed_s000_cer) = iData(tformed_s000_cer) / max(iData(tformed_s000_cer))

DAN_AgLDI_comb <- RegCombIMS::combineReggedIMS(tformed_d001_cer,tformed_s000_cer,mz_offset=2000,sample_name="(-)DAN_AgLDI_comb_2000os")

image(DAN_AgLDI_comb , mz=885.54, main="(-)DAN, combined")
image(DAN_AgLDI_comb , mz=2493.25, main="AgLDI, combined")

#overlay image:
image(DAN_AgLDI_comb , mz=c(885.54,2493.25), col=c("red","green"),
      normalize.image = "linear",superpose=T, main="overlaid ion images", contrast.enhance="suppression")

#See http://cardinalmsi.org for documentation on using Cardinal for data analysis
