# Instance file

This readme contains the explanation of the instance files.

## First Line
The first line has the number of schools, number of students, and maximum number of additional seats as follows:
- nb_schools nb_students nb_max_extra_seats

# School Information (One line per school)
There is one line per school with the name of the school and its initial capacity as follows:
- name_school capacity

# Student-School Pair Information (One line per pair)
There is one line per pair Student-School with the name of the student, the name of the school, the priority score assigned to the student by the school, the utility of this student by attending this school, the probability of acceptance of this student into this school under normal capacity, and the increase in acceptance probability for each additional seat from 1 to nb_max_extra_seats.
- name_student name_school score utility prob_normal_capa inc_prob_install_1_extra_capa inc_prob_install_2_extra_capa ... inc_prob_install_max_extra_capa

# Example Instance
3 5 2  # 3 schools, 5 students, max 2 extra seats
SchoolA 10
SchoolB 8
SchoolC 6
Alice SchoolA 1 8.5 0.6 0.1 0.2
Alice SchoolB 3 7.5 0.5 0.15 0.25
Bob SchoolA 4 8.0 0.55 0.12 0.18
Bob SchoolC 5 7.0 0.4 0.2 0.3
...