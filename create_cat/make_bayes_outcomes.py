#bayes_accept = open('/home/jline/Documents/cataloguing/puma_GLEAM/IDR2/internal_release/puma_gleammulti-v-m-s-n-a-accept.txt').read().split('\n')
#bayes_eyeball = open('/home/jline/Documents/cataloguing/puma_GLEAM/IDR2/internal_release/puma_gleammulti-v-m-s-n-a-accept.txt').read().split('\n')

#bayes_out = open('puma_gleammulti-v-m-s-n-a_outcomes.txt','w+')

#bayes_out.write(bayes_accept[0])

#for line in bayes_accept[1:]:
	#if line != '':
		#bayes_out.write('\n'+line)
	#else:
		#pass
	
#for line in bayes_eyeball:
	#if line != '':
		#bayes_out.write('\n'+line)
	#else:
		#pass
	
#bayes_out.close()

bayes_accept = open('/home/jline/Documents/cataloguing/puma_GLEAM/IDR2/internal_release/puma_gleamdeep-v-m-s-n-a-accept.txt').read().split('\n')
bayes_eyeball = open('/home/jline/Documents/cataloguing/puma_GLEAM/IDR2/internal_release/puma_gleamdeep-v-m-s-n-a-accept.txt').read().split('\n')

bayes_out = open('puma_gleamdeep-v-m-s-n-a_outcomes.txt','w+')

bayes_out.write(bayes_accept[0])

for line in bayes_accept[1:]:
	if line != '':
		bayes_out.write('\n'+line)
	else:
		pass
	
for line in bayes_eyeball:
	if line != '':
		bayes_out.write('\n'+line)
	else:
		pass
	
bayes_out.close()