
class Person():
    def __init__(self, first, last, age, job):
        self.first_name = first
        self.last_name = last
        self.full_name = first+" "+last
        self.age = age
        self.job = job
    def assign_sex(self, sex):
        self.sex = sex


caleb = Person("caleb", "sindic", "27", "programmer")
paul = Person("paul", "runge", "33", "loser")
michelle = Person("michelle","sindic","53","teacher")
people = (caleb, paul, michelle)

for person in people:
    if person.first_name in ("caleb","paul"):
        person.assign_sex("male")
    else:
        person.assign_sex("female")

females = [p.full_name for p in people if p.sex=="female"]
males = [p.full_name for p in people if p.sex=="male"]

print(males,females)
