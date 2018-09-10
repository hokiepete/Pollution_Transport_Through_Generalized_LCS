#!/bin/bash
#This script converts a series of images to a video using ffmpeg
ffmpeg -f image2 -framerate 10 -i 'Vapor_Flux_S1_%04d.tif' -pix_fmt yuv420p -vf scale=1280:-2 'SE_OECS.mp4'
#ffmpeg -f image2 -framerate 10 -i 'SE_lcs_%04d.tif' -pix_fmt yuv420p -vf scale=1280:-2 'SE_LCS.mp4'