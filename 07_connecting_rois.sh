connecting () {

  subj=$1	  

  #original_tracks=/home/auguser2016/Projects/0000_Chiasm/$subj/Tractography_Chiasm_ACT/$subj\_ACT.tck
  #tracks=/home/auguser2016/Projects/0000_Chiasm/$subj/$subj\_Results/$subj\_postLiFE_chiasm.tck
  tracks=$2  

  rois_folder=/home/auguser2016/dMRI_DATA/PREPROCESSED_DATA/Optic_Chiasm_ROIs/$subj
  
  #output_folder=/home/auguser2016/Projects/0000_Chiasm/$subj/$subj\_Results/
  output_folder=$3  

  # load ROIs 
  sr=$rois_folder/$subj\_sr.mif
  sl=$rois_folder/$subj\_sl.mif
  er=$rois_folder/$subj\_er.mif
  el=$rois_folder/$subj\_el.mif

  # Optional: dilate ROIs to make sure that all fibers are classified and perform tracking between dilated ROIs
  
  maskfilter $sr dilate -npass 1 $rois_folder/$subj\_sr_dilated.mif
  maskfilter $sl dilate -npass 1 $rois_folder/$subj\_sl_dilated.mif
  maskfilter $er dilate -npass 1 $rois_folder/$subj\_er_dilated.mif
  maskfilter $el dilate -npass 1 $rois_folder/$subj\_el_dilated.mif
  
  sr_dil=$rois_folder/$subj\_sr_dilated.mif
  sl_dil=$rois_folder/$subj\_sl_dilated.mif
  er_dil=$rois_folder/$subj\_er_dilated.mif
  el_dil=$rois_folder/$subj\_el_dilated.mif
  # perform tracking between ROIs for raw data

  tckedit $tracks -include $sr -include $er $output_folder/$subj\_right_to_right.tck -force
  tckedit $tracks -include $sr -include $el $output_folder/$subj\_right_to_left.tck -force
  tckedit $tracks -include $sl -include $er $output_folder/$subj\_left_to_right.tck -force
  tckedit $tracks -include $sl -include $el $output_folder/$subj\_left_to_left.tck -force
  }

#cd /home/auguser2016/Projects/0000_Chiasm
#for i in */; do
  #connecting $1;
#done

connecting $1 $2 $3
