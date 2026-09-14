

export_selected_audio <- function(wav_file, Samp_rate, bit, highpass) {
      
      if(class(wav_file)!="Wave") {
        wav_file <- Wave(wav_file, samp.rate = Samp_rate, bit = bit)  # Load audio from reactive
      }
      
      
      if (!is.null(wav_file)) {
        
        wav_segment_export <- wav_file

        # Check whether wav_segment is a matrix or a Wave object
        if (inherits(wav_segment_export, "matrix")) {
          # If it's a matrix, we create a Wave object
          wav_segment_export <- Wave(left = as.vector(wav_segment_export), samp.rate = Samp_rate, bit = bit)
        }
        
        # Audio cleaning
        wav_segment_export <- Wave(rmnoise(wav_segment_export, f = Samp_rate), samp.rate = Samp_rate, bit = bit)  # Réduction de bruit
        wav_segment_export <- Wave(fir(wav_segment_export, from = highpass, f = Samp_rate), samp.rate = Samp_rate, bit = bit)  # Eliminer les fréquences en dessous de 250 Hz (High-pass filter)
        
        # normalise to -3db, the standard for eBird bird audio files
        
        # Normalise data by floating [-1,1].
        max_possible <- 2^(bit - 1) - 1
        wav_segment_export@left <- wav_segment_export@left / max_possible
        
        # Calculate the current max
        max_val <- max(abs(wav_segment_export@left))
        
        # Calculate the gain for -3 dB
        target_level <- 10^(-3/20)  # ≈ 0.7079
        gain_norm <- target_level / max_val
        
        # Apply gain
        wav_segment_export@left <- wav_segment_export@left * gain_norm
        
        # Reset data to original scale (integers)
        wav_segment_export@left <- round(wav_segment_export@left * max_possible)
        
        # Convert to mono
        wav_segment_export <- mono(wav_segment_export, which = "left")
        
        wav_segment_export

      } else {
        # Display an error message if it is not a ‘Wave’ object
        showModal(modalDialog(
          title = "Erreur",
          "The audio segment is invalid or not selected.",
          easyClose = TRUE,
          footer = NULL
        ))
      }
    }



