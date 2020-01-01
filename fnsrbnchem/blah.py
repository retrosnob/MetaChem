class Test:
    static_value = 0

    def __init__(self):
        self.value = Test.static_value
        Test.static_value += 1

t1 = Test()
t2 = Test()
for i in range(100):
    t = Test()
print(t.value)

print(t1.value)
print(t2.value)


