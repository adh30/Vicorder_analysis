Augmentation = dataforSV{:,8};
cSBP = dataforSV{:,6};
P1=cSBP-((cSBP-DBP).*Augmentation/100);
ESP = dataforSV{:,15};