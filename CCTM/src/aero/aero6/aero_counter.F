        module aero_counter


        use HDMod
        implicit none

    ! define the counter variables:

        integer, save   :: total_case(3,35), resetted_case(3,35)
        logical         :: flag_criteria1 = .false. , flag_criteria2 = .false.
        real(4)         :: reset_ratio(3,35)

        end module aero_counter