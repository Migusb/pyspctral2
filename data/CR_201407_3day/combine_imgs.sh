
#ffmpeg -r 10 -i ./img/%03d.png -pix_fmt yuv420p -vcodec libx264 out.avi
#ffmpeg -r 10 -i ./img/%03d.png -pix_fmt yuv420p out.avi

convert -delay 30 -loop 1 -dispose Background ./img/*.png combined.gif


