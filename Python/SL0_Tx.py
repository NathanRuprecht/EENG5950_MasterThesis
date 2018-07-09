# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import AudioFile as af
import MicrophoneRecord as mr

def main():
    N=1024
    source="file"
    
    filename = "WhateverItTakes_1min.wav"
    
    timeInterval = 5
    Fs = 8e3
    
    if source=="file":
        af.main(filename, N)
    elif source=="record":
        mr.main(timeInterval, N, Fs)
    else:
        print("No selection")


#Window main parent and title
if __name__ == "__main__":
    main()
