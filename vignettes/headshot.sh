

## https://stackoverflow.com/questions/20075087/how-to-merge-images-in-command-line
# brew install imagemagick 

cd ~/Library/CloudStorage/Dropbox/git/Qindex/vignettes
identify authors.png # 1254x484
identify Amy.png # 214x262 # 214 * 484 / 262 = 395
identify Hallgeir.png # 744x936 # 744 * 484 / 936 = 385
identify Maisel.png # 152x220 # 152 * 484 / 220 = 334
identify Yunguang.png # 212x262 # 212 * 484 / 262 = 392

magick Amy.png -resize 395x484 Amy2.png
magick Hallgeir.png -resize 385x484 Hallgeir2.png
magick Maisel.png -resize 334x484 Maisel2.png
magick Yunguang.png -resize 392x484 Yunguang2.png
magick authors.png Amy2.png Maisel2.png Yunguang2.png Hallgeir2.png +append allauthors.png
