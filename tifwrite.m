function tifwrite(img, savepath)
    fTIF = Fast_Tiff_Write(savepath);
    msg = 0;

    for i = 1:size(img, 3)
        curimg = img(:, :, i)';
        fTIF.WriteIMG(curimg);
    end

    fTIF.close;
end
