Parse_PEAKS_file <- function(path = "~/Desktop/KINATESTID3") {
  
  file_name <- paste(path, "/ABL_PLUS_R1_protein-peptides.csv", sep = "")  
  
  input_file <- read.csv(file_name)
  
  data_rows <- nrow(input_file)
  data_rows <- 100

  #create data frame to capture cases where a single phosphorylation occurred
  single_phos_df <- data.frame(matrix(ncol = 9, nrow = 0))
  
  #create data frame to capture cases where multiple PTMs occurred
  mult_ptms_df <- data.frame(matrix(ncol = 3, nrow = 0))

  # Method 2 to create df  
  # create data frame with 5 empty vectors
  # df2 <- data.frame(Doubles=double(),
  #                   Integers=integer(),
  #                   Factors=factor(),
  #                   Logicals=logical(),
  #                   Characters=character(),
  #                   stringsAsFactors=FALSE)
  
  
  pstr <- input_file$Peptide
  prot_acc <- input_file$Protein.Accession
  
  pep_startl <- as.numeric(input_file$Start)
  seq_startl <- pep_startl 
  
  ascore_str <- input_file$AScore
  ascore_num <- 0
  ascore_split <- ""
  ascore_col <- ""
  central_aa <- ""
  central_aa_desc <- ""
  mod_desc <- ""
  ascore <- ""
  
  i <- 0
  j <- 1
  
  central_aal <- 1
  
  # parensl is short for parenthesis start location
  # parenel is short for parenthsis end locaiton
  parensl = 0
  parenel= 0
  
  # pstr is short for peptide string
  # pstrl is short for peptide string length
  
  # mass_mod represents the ptm in mass of the peptide (in Daltons)
  # mod_type is short for modification, ex: Phosphorylation, Acetylation, etc  
  mass_mod <- ""
  mod_type <- ""
  
  uniprotid <- ""
  protein <- ""
  
  num_ptms <- 0
  num_of_ptms <- 0
  
  while (j <=  data_rows) {   # Beginning of row analysis
    
    pstrl <- nchar(pstr[j])
    parensl <- 0
    parenel <- 0 
    
    # Check to see if this row represents PTM data.  The check is done by 
    # looking for a parenthenthesis in  the string.  The parenthesis will
    # only appear if a PTM did occur.
    # Additionally note: the parenthesis encloses the mass difference in
    # Daltons between a peptide that is phosphorylated and one that is not.
    
    num_of_ptms <- sum(str_count(pstr[j], "\\("))
    
    #######         if (grepl("\\(", pstr[j])) { 
    if (num_of_ptms == 1) {
      
      #      cat("Data Row:")
      #      print(pstr[j])
      
      protinfol <- nchar(prot_acc)  
      #      cat("Protein Accession")
      #      print(prot_acc[j])
      trigger = 0
      
      # This while loop parses out the "UniProt Identifier" and the Protein 
      # from the "Protein Accession" column
      while (i <= protinfol[j]) {
        i = i + 1
        
        if ( "|" == substring(prot_acc[j],i,i) ) {
          trigger = 1
        }
        else {
          if (trigger == 0) {
            uniprotid <- paste(uniprotid, substring(prot_acc[j],i,i), sep = "")
          }
          if (trigger == 1) {
            protein <- paste(protein, substring(prot_acc[j],i,i), sep = "")
          }
        }
      }
      
      i = 0
      
      # This while loop is looking for start and end positions of parentheses
      # and it is also parsing out the mass modificiation number used to 
      # determine the type of ptm that occurred
      while (i <= pstrl) {
        
        i <- i + 1 
        
        if ( "(" == substring(pstr[j],i,i)) {
          parensl = i 
          i = i + 1
          while ( ")" != substring(pstr[j],i,i)) {
            mass_mod <- paste(mass_mod, substring(pstr[j],i,i), sep = "") 
            i = i + 1
          }
          parenel = i  
          i = i + 1
        }
        
      }    # end of while loop looking for left and right parentheses
      
      # if (mass_mod == "+79.97") {
      #   mod_type <- "Phosphoryltation (STY)"
      # }
      # else if (mass_mod == "+.98") {
      #   mod_type <- "Deamidation (NQ)"
      # }
      # else if (mass_mod == "+15.99") {
      #   mod_type <- "Oxidation (M)"
      # }
      # else if (mass_mod == "+17.03") {
      #   mod_type <- "Pyro-glu from Q"
      # }
      # else if (mass_mod == "+57.02") {
      #   mod_type <- "Carbamidomethylation"
      # }
      # else if (mass_mod == "+42.01") {
      #   mod_type <- "Acetylation (Protein N-term)"
      # }
      
      central_aa <- substring(pstr[j], parensl - 1, parensl - 1)
      
      central_aal <- parensl - 1
      
      # Code within this if statement reads the peptide sequence as presented
      # in the output of the PEAKS software program
      
      if (parensl> 0) {    
        
        if (parensl == 4) {
          seq <- "-------"
          seq <- paste(seq, substring(pstr[j],3,3), sep = "") 
        }
        else if (parensl == 5) {
          seq <- "------"
          seq <- paste(seq,substring(pstr[j],3,4),sep="") 
        }
        else if (parensl == 6) {
          seq <- "-----" 
          seq <- paste(seq,substring(pstr[j],3,5),sep="") 
        }
        else if (parensl == 7) {
          seq <- "----" 
          seq <- paste(seq,substring(pstr[j],3,6),sep="") 
        }
        else if (parensl == 8) {
          seq <- "---" 
          seq <- paste(seq,substring(pstr[j],3,7),sep="") 
        }
        else if (parensl == 9) {
          seq <- "--" 
          seq <- paste(seq,substring(pstr[j],3,8),sep="") }
        else if (parensl == 10) {
          seq <- "-" 
          seq <- paste(seq,substring(pstr[j],3,9),sep="") 
        }
        else if (parensl >= 11) {
          seq <- substring(pstr[j], parensl- 8, parensl- 1)  
        }
        
        right_aa_seq_length <- pstrl - parenel- 2
        
        if (right_aa_seq_length == 0)  {
          seq <- paste(seq,"-------",sep="") 
        }
        else if (right_aa_seq_length == 1) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+1), sep="") 
          seq <- paste(seq, "------", sep = "") 
        }
        else if (right_aa_seq_length == 2) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+2), sep="") 
          seq <- paste(seq,"-----",sep="") 
        }
        else if (right_aa_seq_length == 3) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+3), sep="") 
          seq <- paste(seq, "----", sep = "") 
        }
        else if (right_aa_seq_length == 4) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+4), sep="") 
          seq <- paste(seq, "---", sep = "") 
        }
        else if (right_aa_seq_length == 5) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+5), sep="") 
          seq <- paste(seq, "--", sep = "") 
        }
        else if (right_aa_seq_length == 6) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+6), sep="")
          seq <- paste(seq, "-", sep = "") 
        }
        else if (right_aa_seq_length  >= 7) {
          seq <- paste(seq, substring(pstr[j], parenel+1, parenel+7), sep="") 
        }
        
      }
      
      # cat("Row Number")
      # print(j)
      # cat("File Data Length") 
      # print(pstrl)    
      # cat("File Parenthesis Start") 
      # print(parensl)
      # cat("File Parenthesis End") 
      # print(parenel)
      # cat("Central Amino Acid") 
      # print(central_aa)
      # cat("Mass Modification") 
      # print(mass_mod)
      # cat("Modification Type") 
      # print(mod_type)
      #      cat("Amino Acid Sequence") 
      #      print(seq)
      # print("UniProt ID")
      # cat("Uniprotid")
      # print(uniprotid)
      # cat("Protein")
      # print(protein)     
      # cat("\r\n")
      
      seq_startl = pep_startl[j] + central_aal - 3
      
      ascore_num <- str_count(ascore_str[j], ";")
      if (ascore_num == 0) {
        ascore_split <- ascore_str[j]
      }
      else {
        ascore_split <- str_split(ascore_str[j], ";", simplify = TRUE)    
      }
      # cat("AScore String\r\n")
      k = 1
      while (k <= (ascore_num+1)) {
        ascore_col <- str_split(ascore_split[1], ":", simplify = TRUE)
        mod_desc <- ascore_col[2]
        ascore <- ascore_col[3]
        # print(j)
        # print(uniprotid)
        # print(protein)
        # print(central_aa)
        # print(mass_mod)
        # print(mod_desc)
        # print(ascore)     
        # print(seq_startl)
        # print(seq)
        k = k + 1
      }
      # cat("Central AA\r\n")
      #central_aa_desc <- ascore_col[1]
      #          print(central_aa_desc)
      # cat("Modification Description\r\n")
      #         mod_desc <- ascore_col[2]
      #         print(mod_desc)
      #          cat("AScore\r\n")
      
      #print(mod_type)
      
      data_vector <- c(j,uniprotid,protein,central_aa,mass_mod,mod_desc,ascore,     
                       seq_startl,seq)
      single_phos_df <- rbind(single_phos_df, data_vector)
      # print(data_vector)
      
      # cat("\r\n")
    }   # End of sum(grep) == 1 to confirm 1 PTM in the row
    else if (num_of_ptms > 1) {
#      cat("MULTIPLE PTMS in this ROW\r\n")
      mult_ptms_v <- c(pstr[j],j,num_of_ptms)
      mult_ptms_df <- rbind(mult_ptms_df,mult_ptms_v)
      # print(mult_ptms_df)
    }
    
    i <- 0
    j <- j + 1  
    seq <- ""
    mod_type <- ""
    mass_mod <- ""
    parensl = 0
    parenel= 0
    uniprotid <- ""
    protein <- ""
    ascore_num <- 0
    ascore_split <- ""
    ascore_col <- ""
    central_aa_desc <- ""
    mod_desc <- ""
    ascore <- ""
    num_of_ptms <- 0
    
    
    
  }     # End of row analysis
  
  colnames(mult_ptms_df) <- c("Multiple PTMs","Row Number", "Number of PTMs")
  colnames(single_phos_df) <- c("Row", "UniProt Identifier", "Protein", 
                                "Amino Acid Letter", "Modification Number", 
                                "Modification Name", "AScore", "Position",
                                "Peptide")
  print(mult_ptms_df)
  print(single_phos_df)
  
}   # End of function
