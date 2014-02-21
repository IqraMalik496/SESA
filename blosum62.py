#CODE CHALLENGE: Solve the Global Alignment Problem.
#     Input: Two protein strings written in the single-letter amino acid alphabet.
#     Output: The maximum alignment score of these strings followed by an alignment achieving this
#     maximum score. Use the BLOSUM62 scoring matrix and indel penalty sigma = 5.

Input="""
PLEASANTLY
MEANLY
"""
Input="""
KMMDNKEDLKAVEWHIQHGRCKQTQALEHGTLAEVPHHAWLQNPWDIAISSDASVIKWIWPDMTHNCVHEGPDDFRFVMCPGYMTHSVLPQVHRPVLWIKGSPGDLSWYQWTIHWYTHCMVNHGSCRKVMDWRTWCVYMILECWAFMVCTFQFWYPNTRCHCDWAMMRLCFDGFLRKMNGAWSVIRCPQYKCVWTTWFLQKGNFCRYWMVWVFVCWNELISYFTTMIKFWDDNELFSGYGCWIWFENGTSWNLWLMTKHFQRQCNSSIMRESWYEGPDNSVIGRYGYKIIDKFSPLHPQSPYMEFRENTLVCMIYQLYHHNCGPRARKPRLWAFTVGVMCLPIMKLIPLKDKRKQVSAEDHLNFSRDAEITRRFEIGWLCFAKHHMIGTRWFPADMKMGWMYCEQTLTWQLLYPVHLKLWRGLKSSISAVRPGYGKYRTQNSKYSYCMMCHQQVANGITINRGNYMLFLSCVKEMDARSRVETCSEHRMQVKQTIREKIGKDYQPSGKCCKFCCTMHRIFCSVENRGATYYVAEMPVPYSSNDAFHDSWVDFWFAEQGVMQKSYHPVKWKYHVDLMGLIDSFCHDVMKANPLLEIGYRAHCHSHYAQNRVHYMPFSYPDTDVIEWISGYVENSGLRKYIENQTYVHGWSNYPKTLCNSKNEEPIPGHVKWCAMKCGVAVSQYRSVQMYEVVENIHSKWERKCYRWHKFAIQLALMTQRENMDNWDVCFLMFNQRFGTSHPEAQFRYITWHRDCRLSQFDCKQLCARNCPTRQGKTPGQVRKVITY
KQIMDMDYKNKEDCYRATLKLESMAVEDHMVGIQHQTTQRCNWYKFQTQASEWWEHGTLAEVPHEKDSKLDMWLQNPWDIAHSSDADMTHNCVHEGPDDFRFVMWPNYMTHSVVLWIKGWTIHWYTHCNYNHGSCRKVMDWRTWCVYMILVVEHGWTFQFEFVRWPNTRCHADWAMMRLCFDGFLYKLGAWSVHEHPFKNMAKCVTKGNFCRYWMVIVFVHVTSIFLISTMTTMIKNWDDNELFSGYGCCQPFGIWFENGDSWNLWLMTKHFQRQCKSSIMRESWYATTWSDSVIGRYGILHPPFWHENTCVCIYQLRARKAGEIKHVVMCLPIMKIPLKSKRKQVSREDHLNFNTLDAWRDQEMWFDNEYLNTRRFEIGWLCMFKHHLIGTRWFPADIKMNKMGWMERCLPVHVLSEAQSRKQWRMQRSCGYGRSVWNYYTTYGTQNSKYSYCMVVCHQQVATICRGNKEMDARSRVATCSEHRETKFLMHGKMVCTPKPCSGKCCKFCCTMHRIFCSVENRGATYYVAEMPVPYSNIAFHMNSICQFWWARMDMWNNHLHVALMGLMKANPLLEIGYRAHCSHYAGYRGNRVHYMPFSYPDTDCIECISGYHENSGWMKYIENQTPCVHGWSNCNAKNEEPIPGHVKKCGVAPYVEGCSASRSVQMYEVVQKCYRWHKMMGSAIQGRYIVMARMTQRENMDEWDVCFLMFRFGTSHPEAQFRYILDCRLSQFDPKQLCARNCPTRQGKTPGQVRKVITY
"""
#
#Sample Output:
#     8
#     PLEASANTLY
#     -MEA--N-LY

blosum62="""
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2
V  0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
"""
blosum62= [line for line in blosum62.splitlines() if line]
Amino_acids = blosum62[0].split()
blosum_dict = {}
for i in range(len(Amino_acids)):
    for j in range(len(Amino_acids)):
        blosum_dict[Amino_acids[i]+Amino_acids[j]]=blosum62[i+1].split()[j+1]
#print blosum_dict

[v,w]=Input.split()

def outputLCS(backtrack,v,i,j,v_list,w_list):
    if i==0 and j==0:
        return
    if backtrack[i,j]=='south':
        outputLCS(backtrack,v,i-1,j,v_list,w_list)
        v_list.append(v[i-1])
        w_list.append('-')
    elif backtrack[i,j]=='east':
        outputLCS(backtrack,v,i,j-1,v_list,w_list)
        v_list.append('-')
        w_list.append(w[j-1])
    else:
        outputLCS(backtrack,v,i-1,j-1,v_list,w_list)
#        print v[i]
        v_list.append(v[i-1])
        w_list.append(w[j-1])

def lcs(v,w):
    sigma=5
    s={}
    backtrack={}
    for i in range(len(v)+1):
        s[i,0]=-i*sigma
        backtrack[i,0]='south'
    for j in range(1,len(w)+1):
        s[0,j]=-j*sigma
        backtrack[0,j]='east'
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            s[i,j]=max(s[i-1,j]-sigma,
                       s[i,j-1]-sigma,
                       s[i-1,j-1]+int(blosum_dict[v[i-1]+w[j-1]]))
            backtrack[i,j]='south' if (s[i,j]==s[i-1,j]-sigma) else 'east' if (s[i,j]==s[i,j-1]-sigma) else 'southeast'
    return s, backtrack

if __name__=='__main__':
    import sys
    sys.setrecursionlimit(10000)
    scores,backtrack=lcs(v,w)
    v_list=[]
    w_list=[]
    outputLCS(backtrack,v,len(v),len(w),v_list,w_list)
#print v
#print w

#print scores[len(v),len(w)]
#print ''.join(v_list)
#print ''.join(w_list)
#for i in range(len(v)+1):
#    print [scores[i,j] for j in range(len(w)+1)]
