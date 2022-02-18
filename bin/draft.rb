def divide_pair_into_n(globListData, nAmount, fwdIdentifier, revIdentifier)
  allData=Dir.glob(globListData)
  pair_ch=[]
  allData.each do |x|
    if x.match("_#{fwdIdentifier}")
      prefix=x.match("_#{fwdIdentifier}").pre_match
      suffix=x.match("_#{fwdIdentifier}").post_match
      pair_ch.push("#{prefix}{#{fwdIdentifier},#{revIdentifier}}#{suffix}")
    end
  end
  equalarray=pair_ch.each_slice(nAmount).to_a.count

  chunkArray=pair_ch.each_slice(equalarray).to_a
  for i in 0..chunkArray.count-1
    write=File.open("#{File.dirname(__FILE__)}/list#{i+1}.txt", "w")
    for k in 0..chunkArray[i].count-1
      write.puts(chunkArray[i][k])
    end
  end
end
