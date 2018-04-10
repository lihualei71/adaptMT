for (file in list.files()){
    if (file != "source-all.R"){
        source(file)
    }
}
