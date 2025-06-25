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
  
## Description
# Hemoglobin concentration at baseline, 4 weeks and 8 weeks in carcinomatous
# anemia patients from a randomized clinical trial.

## Variables
# X_1: Baseline, 
# X_2: 4 weeks, 
# X_3: 8 weeks

## Categories
# ≥ 10 g/dl
# 8−10 g/dl 
#  < 8 g/dl

database$hemoglobin <- array(
  c(
    77, 43, 3,
    3,  17, 3,
    1,   0, 0,
    
    7,  7, 0,
    8, 16, 8,
    1,  2, 4,
    
    1, 0, 0,
    1, 5, 1,
    1, 3, 3
  ),
  dim = c(3, 3, 3)
)


## Description
# Stationary two-step transitions in a panel study of potential voters
# in Erie County, Ohio, 1940 (from Bishop et al., 1975, p.305).

## Variables
# X_1: Time t-2, 
# X_2: Time t-1, 
# X_3: Time  t

## Categories
# Republican,
# Democrat,
# Undecided

database$voters <- array(
  c(
    557, 71, 18,
    17,  62,  4,
    3,    6,  9, 
    
    16,  11,  5,
    21, 346, 24,
    0,    6, 22,
    
    6,   1,   0,
    5,  54,  10,
    8,  63, 435 
  ),
  dim = c(3, 3, 3)
)


# Description
# Crossover study for treating dysmenorrhea (Kenward and Jones, 1991).

## Variables
# X_1: Placebo,
# X_2: Low-dose analgesic, 
# X_3: High-dose analgesic

## Categories
# giving no relief,
# moderate relief,
# complete relief

database$dysmenorrhea <- array(
  c(
    6, 2, 1, 
    3, 1, 0,
    1, 2, 1,
    
    4, 3, 0, 
    13, 3, 0,
    8, 1, 1,
    
    5, 2, 2,
    10, 1, 0,
    14, 2, 0
  ),
  dim = c(3, 3, 3)
)


# Description
# Opinions about government spending in 2022 from the GSS (Yoshimoto et al., 2020).

## Variables
# X_1: Education,
# X_2: Environment, 
# X_3: Assistance to the poor

## Categories
# too little,
# about right,
# too much 

database$opinions1 <- array(
  c(
    612, 85, 12,
    134, 46, 16,
     51,  9, 13,
    
    110, 30,  8,
     55, 43, 16,
     11, 11,  8,
    
    30, 6, 3,
    11, 9, 8,
    11, 5, 13
  ),
  dim = c(3, 3, 3)
)


## Description
# Opinions about government spending in 2022 from the GSS.

## Variables
# X_1: National defense, 
# X_2: Space exploration, 
# X_3: Health

## Categories
# too little,
# about right,
# too much 

database$opinions2 <- array(
  c(
    88,  85,  77,
    178, 261, 127,
    108,  89, 110,
    
    24,  14, 15,
    50, 114, 23,
    14,  19, 18,
    
    33, 10,  7, 
    20, 48, 19,
    11, 25,  8
  ),
  dim = c(3, 3, 3)
)

## Description
# Which comes closest to your view about what government policy should be toward
# unauthorized immigrants now living in the United States?

## Variables
# X_1: wave 1
# X_2: wave 2
# X_3: wave 3

## Categories
# (1) Make all unauthorized immigrants felons and send them back to their home country.
# (2) Have a guest worker program that allows unauthorized immigrants to remain in
#     the United States in order to work, but only for a limited time.
# (3) Allow unauthorized immigrants to remain in the United States and eventually
#     qualify for U.S. citizenship, but only if they meet certain conditions.

database$immigrants <- array(
  c(
    300,  61,  34,
     66,  47,  25,
     29,  20,  18,
    
     68,  49,  21,
    108, 324, 136,
     36, 135, 166,
    
     35,  22,   33,
     42, 123,  173,
     88, 342, 1775
  ),
  dim = c(3, 3, 3)
)


## Description
# Do you think of yourself as closer to the Republican Party or to the Democratic Party?

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