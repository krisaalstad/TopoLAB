clearvars;
urlis='https://www.dropbox.com/s/pmp5nobqj9cmvp1/input.zip?dl=1'; % Changed dl=0 to dl=1
tarf='input.zip';
websave(tarf,urlis);
unzip(tarf,'tmp') 