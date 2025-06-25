#
# Sourcing this R file contains database used in main.R.
#
# Create an environment to store all datasets
  database <- new.env(parent = emptyenv())
#
#
# database$example <- array(
#  c(
#    (1,1,1), (2,1,1), (3,1,1),
#    (1,2,1), (2,2,1), (3,2,1),
#    (1,3,1), (2,3,1), (3,3,1),
#
#    (1,1,2), (2,1,2), (3,1,2),
#    (1,2,2), (2,2,2), (3,2,2),
#    (1,3,2), (2,3,2), (3,3,2),
#    
#    (1,1,3), (2,1,3), (3,1,3),
#    (1,2,3), (2,2,3), (3,2,3),
#    (1,3,3), (2,3,3), (3,3,3),
#  ),
#  dim = c(3, 3, 3)
#  )
#
#########################  BEGIN import database  ##############################

## Description
# Political identification across three waves of the ANES 2020–2022 Social Media Study.
# Questionnaire: "Do you think of yourself as closer to the Republican Party or to the Democratic Party?"

## Variables
# X_1: wave 1
# X_2: wave 2
# X_3: wave 3

## Categories
# (1) Closer to the Republican Party
# (2) Neither
# (3) Closer to the Democratic Party

database$party <- array(
  c(
    240, 20,   4,
    11, 18,   0,
    0,  1,   7,
    
    32,  22,   0,
    23, 237,  28,
    2,  24,  36,
    
    8,  4,   5,
    5, 28,  16, 
    4, 29, 323
  ),
  dim = c(3, 3, 3)
)
  
##########################  END import database  ###############################
